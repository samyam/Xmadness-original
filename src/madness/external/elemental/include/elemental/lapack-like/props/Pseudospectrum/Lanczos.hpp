/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_LANCZOS_HPP
#define ELEM_PSEUDOSPECTRUM_LANCZOS_HPP

#include "./Power.hpp"

namespace elem {
namespace pspec {

const Int HCapacityInit = 10;

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real>>& HDiagList, 
  const std::vector<std::vector<Real>>& HSubdiagList,
  Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;
    const Int krylovSize = HDiagList[0].size();
    std::vector<Real> HDiag, HSubdiag, w(krylovSize);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j]; 
        HSubdiag = HSubdiagList[j];
        if( !HasNan(HDiag) && !HasNan(HSubdiag) )
        {
            lapack::SymmetricTridiagEig     
            ( krylovSize, HDiag.data(), HSubdiag.data(), w.data(), 
              krylovSize-1, krylovSize-1 );
            const Real est = Sqrt(w[0]);
            activeEsts.Set( j, 0, Min(est,normCap) );
        }
        else
            activeEsts.Set( j, 0, normCap );
    }
}

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real>>& HDiagList, 
  const std::vector<std::vector<Real>>& HSubdiagList,
  DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts.Matrix() );
}

template<typename Real>
inline void
Deflate
( std::vector<std::vector<Real>>& HDiagList,
  std::vector<std::vector<Real>>& HSubdiagList,
  Matrix<Complex<Real>>& activeShifts, 
  Matrix<Int          >& activePreimage,
  Matrix<Complex<Real>>& activeXOld,
  Matrix<Complex<Real>>& activeX,
  Matrix<Real         >& activeEsts, 
  Matrix<Int          >& activeConverged,
  Matrix<Int          >& activeItCounts,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    Timer timer;
    if( progress )
        timer.Start();
    const Int numActive = activeX.Width(); 
    Int swapTo = numActive-1;
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( activeConverged.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                std::swap( HDiagList[swapFrom], HDiagList[swapTo] );
                std::swap( HSubdiagList[swapFrom], HSubdiagList[swapTo] );
                RowSwap( activeShifts,   swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts,     swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                ColSwap( activeXOld,     swapFrom, swapTo );
                ColSwap( activeX,        swapFrom, swapTo );
            }
            --swapTo;
        }
    }
    if( progress )
        std::cout << "Deflation took " << timer.Stop() << " seconds"
                  << std::endl;
}

template<typename Real>
inline void
Deflate
( std::vector<std::vector<Real>>& HDiagList,
  std::vector<std::vector<Real>>& HSubdiagList,
  DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Complex<Real>        >& activeXOld,
  DistMatrix<Complex<Real>        >& activeX,
  DistMatrix<Real,         MR,STAR>& activeEsts,
  DistMatrix<Int,          MR,STAR>& activeConverged,
  DistMatrix<Int,          VR,STAR>& activeItCounts,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    Timer timer;
    if( progress && activeShifts.Grid().Rank() == 0 )
        timer.Start();
    const Int numActive = activeX.Width(); 
    Int swapTo = numActive-1;

    DistMatrix<Complex<Real>,STAR,STAR> shiftsCopy( activeShifts );
    DistMatrix<Int,STAR,STAR> preimageCopy( activePreimage );
    DistMatrix<Real,STAR,STAR> estimatesCopy( activeEsts );
    DistMatrix<Int, STAR,STAR> itCountsCopy( activeItCounts );
    DistMatrix<Int, STAR,STAR> convergedCopy( activeConverged );
    DistMatrix<Complex<Real>,VC,STAR> XOldCopy( activeXOld ), XCopy( activeX );

    const Int n = ( activeX.LocalWidth()>0 ? HDiagList[0].size() : 0 );
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedCopy.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                // TODO: Avoid this large latency penalty
                if( activeX.IsLocalCol(swapFrom) && 
                    activeX.IsLocalCol(swapTo) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    const Int localTo = activeX.LocalCol(swapTo);
                    DEBUG_ONLY(
                        if( HDiagList[localFrom].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HDiagList[localTo].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HSubdiagList[localFrom].size() != n )
                            LogicError("Invalid HSubdiagList size");
                        if( HSubdiagList[localTo].size() != n )
                            LogicError("Invalid HSubdiagList size");
                    )
                    std::swap( HDiagList[localFrom], HDiagList[localTo] );
                    std::swap( HSubdiagList[localFrom], HSubdiagList[localTo] );
                }
                else if( activeX.IsLocalCol(swapFrom) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    DEBUG_ONLY(
                        if( HDiagList[localFrom].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HSubdiagList[localFrom].size() != n )
                            LogicError("Invalid HSubdiagList size");
                    )
                    const Int partner = activeX.ColOwner(swapTo);
                    mpi::TaggedSendRecv
                    ( HDiagList[localFrom].data(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiagList[localFrom].data(), n, 
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }
                else if( activeX.IsLocalCol(swapTo) )
                {
                    const Int localTo = activeX.LocalCol(swapTo);
                    DEBUG_ONLY(
                        if( HDiagList[localTo].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HSubdiagList[localTo].size() != n )
                            LogicError("Invalid HSubdiagList size");
                    )
                    const Int partner = activeX.ColOwner(swapFrom);
                    mpi::TaggedSendRecv
                    ( HDiagList[localTo].data(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiagList[localTo].data(), n, 
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }

                RowSwap( shiftsCopy,    swapFrom, swapTo );
                RowSwap( preimageCopy,  swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy,  swapFrom, swapTo );
                ColSwap( XOldCopy,      swapFrom, swapTo );
                ColSwap( XCopy,         swapFrom, swapTo );
            }
            --swapTo;
        }
    }

    activeShifts   = shiftsCopy;
    activePreimage = preimageCopy;
    activeEsts     = estimatesCopy;
    activeItCounts = itCountsCopy;
    activeXOld     = XOldCopy;
    activeX        = XCopy;

    if( progress ) 
    {
        mpi::Barrier( activeShifts.Grid().Comm() );
        if( activeShifts.Grid().Rank() == 0 ) 
            std::cout << "Deflation took " << timer.Stop() << " seconds"
                      << std::endl;
    }
}

template<typename Real>
inline Matrix<Int>
Lanczos
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Lanczos"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();

    const Int maxIts = psCtrl.maxIts;
    const bool deflate = psCtrl.deflate;
    const bool progress = psCtrl.progress;

    // Keep track of the number of iterations per shift
    Matrix<Int> itCounts;
    Ones( itCounts, numShifts, 1 );

    // Keep track of the pivoting history if deflation is requested
    Matrix<Int> preimage;
    Matrix<C> pivShifts( shifts );
    if( deflate )
    {
        preimage.Resize( numShifts, 1 );
        for( Int j=0; j<numShifts; ++j )
            preimage.Set( j, 0, j );
    }

    // The Hessenberg case currently requires explicit access to the adjoint
    Matrix<C> UAdj, activeShiftsConj;
    if( !psCtrl.schur )
        Adjoint( U, UAdj );

    // Simultaneously run Lanczos for various shifts
    Matrix<C> XOld, X, XNew;
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiagList( numShifts ),
                                   HSubdiagList( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        HDiagList[j].reserve( HCapacityInit );
        HSubdiagList[j].reserve( HCapacityInit-1 );
    }

    psCtrl.snapCtrl.ResetCounts();

    Timer timer, subtimer;
    Int numIts=0, numDone=0;
    Matrix<Real> estimates(numShifts,1);
    Zeros( estimates, numShifts, 1 );
    auto lastActiveEsts = estimates;
    Matrix<Int> activePreimage;
    std::vector<Real> realComponents;
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = View( pivShifts, 0, 0, numActive, 1 );
        auto activeEsts = View( estimates, 0, 0, numActive, 1 );
        auto activeItCounts = View( itCounts, 0, 0, numActive, 1 );
        auto activeXOld = View( XOld, 0, 0, n, numActive );
        auto activeX    = View( X,    0, 0, n, numActive );
        auto activeXNew = View( XNew, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );
        HDiagList.resize( numActive );
        HSubdiagList.resize( numActive );

        if( progress )
            timer.Start();
        activeXNew = activeX;
        if( psCtrl.schur )
        {
            if( progress )
                subtimer.Start();
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeXNew );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (8.*n*n*numActiveShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftTrsm's: " << msTime << " seconds, "
                          << gflops << " GFlops" << std::endl;
            }
        }
        else
        {
            if( progress )
                subtimer.Start();
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj, activeShiftsConj, activeXNew );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (32.*n*n*numActiveShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftHessSolve's: " << msTime 
                          << " seconds, " << gflops << " GFlops" << std::endl;
            }
        }

        // Orthogonalize with respect to the old iterate
        if( numIts > 0 )
        {
            ExtractList( HSubdiagList, realComponents, numIts-1 );
            ColumnSubtractions( realComponents, activeXOld, activeXNew );
        }

        // Orthogonalize with respect to the last iterate
        InnerProducts( activeX, activeXNew, realComponents );
        PushBackList( HDiagList, realComponents );
        ColumnSubtractions( realComponents, activeX, activeXNew );

        // Compute the norm of what is left
        ColumnNorms( activeXNew, realComponents );
        PushBackList( HSubdiagList, realComponents );

        activeXOld = activeX;
        activeX    = activeXNew; 
        InvBetaScale( realComponents, activeX );
        if( progress )
            subtimer.Start();
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
            std::cout << "  Ritz computations: " << subtimer.Stop() 
                      << " seconds" << std::endl;

        auto activeConverged = 
            FindConverged
            ( lastActiveEsts, activeEsts, activeItCounts, psCtrl.tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress )
        {
            const double iterTime = timer.Stop();
            std::cout << "iteration " << numIts << ": " << iterTime
                      << " seconds, " << numDone << " of " << numShifts
                      << " converged" << std::endl;
        }

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
            Deflate
            ( HDiagList, HSubdiagList, activeShifts, activePreimage, activeXOld,
              activeX, activeEsts, activeConverged, activeItCounts, progress );

        lastActiveEsts = activeEsts;

        // Save snapshots of the estimates at the requested rate
        psCtrl.snapCtrl.Iterate();
        Snapshot
        ( preimage, estimates, itCounts, numIts, deflate, psCtrl.snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );
    FinalSnapshot( invNorms, itCounts, psCtrl.snapCtrl );

    return itCounts;
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
Lanczos
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Lanczos"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = U.Grid();

    const Int maxIts = psCtrl.maxIts;
    const bool deflate = psCtrl.deflate;
    const bool progress = psCtrl.progress;

    if( deflate && g.Rank() == 0 ) 
        std::cerr << "NOTE: Deflation swaps not yet optimized!" << std::endl;

    // Keep track of the number of iterations per shift
    DistMatrix<Int,VR,STAR> itCounts(g);
    Ones( itCounts, numShifts, 1 );

    // Keep track of the pivoting history if deflation is requested
    DistMatrix<Int,VR,STAR> preimage(g);
    DistMatrix<C,  VR,STAR> pivShifts( shifts );
    if( deflate )
    {
        preimage.AlignWith( shifts );
        preimage.Resize( numShifts, 1 );
        const Int numLocShifts = preimage.LocalHeight();
        for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
        {
            const Int i = preimage.GlobalRow(iLoc);
            preimage.SetLocal( iLoc, 0, i );
        }
    }

    // The Hessenberg case currently requires explicit access to the adjoint
    DistMatrix<C,VC,STAR> U_VC_STAR(g), UAdj_VC_STAR(g);
    DistMatrix<C,VR,STAR> activeShiftsConj(g);
    DistMatrix<C,STAR,VR> activeXNew_STAR_VR(g);
    if( !psCtrl.schur )
    {
        U_VC_STAR = U;
        Adjoint( U, UAdj_VC_STAR );
    }

    // Simultaneously run Lanczos for various shifts
    DistMatrix<C> XOld(g), X(g), XNew(g);
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiagList( X.LocalWidth() ),
                                   HSubdiagList( X.LocalWidth() );
    for( Int j=0; j<X.LocalWidth(); ++j )
    {
        HDiagList[j].reserve( HCapacityInit );
        HSubdiagList[j].reserve( HCapacityInit-1 );
    }

    psCtrl.snapCtrl.ResetCounts();

    Timer timer, subtimer;
    Int numIts=0, numDone=0;
    DistMatrix<Real,MR,STAR> estimates(g);
    estimates.AlignWith( shifts );
    Zeros( estimates, numShifts, 1 );
    auto lastActiveEsts = estimates;
    DistMatrix<Int,VR,STAR> activePreimage(g);
    std::vector<Real> realComponents;
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = View( pivShifts, 0, 0, numActive, 1 );
        auto activeEsts = View( estimates, 0, 0, numActive, 1 );
        auto activeItCounts = View( itCounts, 0, 0, numActive, 1 );
        auto activeXOld = View( XOld, 0, 0, n, numActive );
        auto activeX    = View( X,    0, 0, n, numActive );
        auto activeXNew = View( XNew, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );
        HDiagList.resize( activeX.LocalWidth() );
        HSubdiagList.resize( activeX.LocalWidth() );

        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                timer.Start();
        }
        activeXNew = activeX;
        if( psCtrl.schur )
        {
            if( progress )
            { 
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                    subtimer.Start();
            }
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeXNew );
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = (8.*n*n*numActiveShifts)/(msTime*1e9);
                    std::cout << "  MultiShiftTrsm's: " << msTime 
                              << " seconds, " << gflops << " GFlops" 
                              << std::endl;
                }
            }
        }
        else
        {
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                    subtimer.Start();
            }
            // NOTE: This redistribution sequence might not be necessary
            activeXNew_STAR_VR = activeXNew;
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U_VC_STAR, activeShifts,
              activeXNew_STAR_VR );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj_VC_STAR, activeShiftsConj,
              activeXNew_STAR_VR );
            activeXNew = activeXNew_STAR_VR;
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = 
                        (32.*n*n*numActiveShifts)/(msTime*1.e9);
                    std::cout << "  MultiShiftHessSolve's: " << msTime
                              << " seconds, " << gflops << " GFlops" 
                              << std::endl;
                }
            }
        }

        // Orthogonalize with respect to the old iterate
        if( numIts > 0 )
        {
            ExtractList( HSubdiagList, realComponents, numIts-1 );
            ColumnSubtractions( realComponents, activeXOld, activeXNew );
        }

        // Orthogonalize with respect to the last iterate
        InnerProducts( activeX, activeXNew, realComponents );
        PushBackList( HDiagList, realComponents );
        ColumnSubtractions( realComponents, activeX, activeXNew );

        // Compute the norm of what is left
        ColumnNorms( activeXNew, realComponents );
        PushBackList( HSubdiagList, realComponents );

        activeXOld = activeX;
        activeX    = activeXNew;
        InvBetaScale( realComponents, activeX );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                subtimer.Start();
        }
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                std::cout << "  Ritz computations: " << subtimer.Stop() 
                          << " seconds" << std::endl;
        }

        auto activeConverged =
            FindConverged
            ( lastActiveEsts, activeEsts, activeItCounts, psCtrl.tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
            {
                const double iterTime = timer.Stop();
                std::cout << "iteration " << numIts << ": " << iterTime
                          << " seconds, " << numDone << " of " << numShifts
                          << " converged" << std::endl;
            }
        }

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
            Deflate
            ( HDiagList, HSubdiagList, activeShifts, activePreimage, activeXOld,
              activeX, activeEsts, activeConverged, activeItCounts, progress );

        lastActiveEsts = activeEsts;

        // Save snapshots of the estimates at the requested rate
        psCtrl.snapCtrl.Iterate();
        Snapshot
        ( preimage, estimates, itCounts, numIts, deflate, psCtrl.snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );
    FinalSnapshot( invNorms, itCounts, psCtrl.snapCtrl );

    return itCounts;
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_LANCZOS_HPP
