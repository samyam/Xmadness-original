/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_BUSINGERGOLUB_HPP
#define ELEM_QR_BUSINGERGOLUB_HPP

#include ELEM_DIAGONALSCALETRAPEZOID_INC
#include ELEM_GEMV_INC
#include ELEM_GER_INC
#include ELEM_SWAP_INC

#include ELEM_REFLECTOR_INC

#include ELEM_INVERTPERMUTATION_INC

#include ELEM_ZEROS_INC

#include <algorithm>

namespace elem {
namespace qr {

template<typename F>
inline Base<F>
ColNorms( const Matrix<F>& A, std::vector<Base<F>>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ColNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real maxNorm = 0;
    norms.resize( n );
    for( Int j=0; j<n; ++j )
    {
        norms[j] = blas::Nrm2( m, A.LockedBuffer(0,j), 1 );
        maxNorm = std::max( maxNorm, norms[j] );
    }
    return maxNorm;
}

template<typename Real>
inline ValueInt<Real>
FindPivot( const std::vector<Real>& norms, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("qr::FindPivot"))
    const auto maxNorm = std::max_element( norms.begin()+col, norms.end() );
    ValueInt<Real> pivot;
    pivot.value = *maxNorm;
    pivot.index = maxNorm - norms.begin();
    return pivot;
}

template<typename F> 
inline Int
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, Matrix<Int>& pPerm,
  Int maxSteps, Base<F> tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    maxSteps = Min(maxSteps,Min(m,n));
    t.Resize( maxSteps, 1 );
    d.Resize( maxSteps, 1 );

    Matrix<F> z21;

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms;
    const Real maxOrigNorm = ColNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());

    // Initialize the inverse permutation to the identity
    Matrix<Int> pInvPerm;
    pInvPerm.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
        pInvPerm.Set( j, 0, j ); 

    Int k=0;
    for( ; k<maxSteps; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find the next column pivot
        const ValueInt<Real> pivot = FindPivot( norms, k );
        if( pivot.value <= tol*maxOrigNorm )
            break;
        RowSwap( pInvPerm, k, pivot.index );
 
        // Perform the swap
        const Int jPiv = pivot.index;
        if( jPiv != k )
        {
            blas::Swap( m, A.Buffer(0,k), 1, A.Buffer(0,jPiv), 1 );
            norms[jPiv] = norms[k];
            origNorms[jPiv] = origNorms[k];
        }

        // Find tau and u such that
        //  / I - tau | 1 | | 1, u^H | \ | alpha11 | = | beta |
        //  \         | u |            / |     a21 | = |    0 |
        const F tau = LeftReflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        // AB2 := Hous(aB1,tau) AB2
        //      = (I - tau aB1 aB1^H) AB2
        //      = AB2 - tau aB1 (AB2^H aB1)^H
        Zeros( z21, AB2.Width(), 1 );
        Gemv( ADJOINT, F(1), AB2, aB1, F(0), z21 );
        Ger( -tau, aB1, z21, AB2 );

        // Reset alpha11's value
        alpha11.Set(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176
        for( Int j=k+1; j<n; ++j )
        {
            if( norms[j] != Real(0) )
            {
                Real gamma = Abs(A.Get(k,j)) / norms[j];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );

                const Real ratio = norms[j] / origNorms[j];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || alwaysRecompute )
                {
                    norms[j] = blas::Nrm2( m-(k+1), A.Buffer(k+1,j), 1 );
                    origNorms[j] = norms[j];
                }
                else
                    norms[j] *= Sqrt(gamma);
            }
        }
    }
    InvertPermutation( pInvPerm, pPerm );

    // Form d and rescale R
    auto R = View( A, 0, 0, k, n );
    d = R.GetRealPartOfDiagonal();
    for( Int j=0; j<k; ++j )
    {
        const Real delta = d.Get(j,0);
        if( delta >= Real(0) )
            d.Set(j,0,Real(1));
        else
            d.Set(j,0,Real(-1));
    }
    DiagonalScaleTrapezoid( LEFT, UPPER, NORMAL, d, R );

    return k;
}

template<typename F> 
inline Int
BusingerGolub
( Matrix<F>& A, Matrix<Int>& pPerm,
  Int maxSteps, Base<F> tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    return BusingerGolub( A, t, d, pPerm, maxSteps, tol, alwaysRecompute );
}

template<typename F> 
inline Int
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, Matrix<Int>& pPerm,
  Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    // Use a tolerance of -1 so that we do not stop early
    return BusingerGolub
           ( A, t, d, pPerm, numSteps, Base<F>(-1), alwaysRecompute );
}

// If we don't need 't' or 'd' from the above routine
template<typename F> 
inline Int
BusingerGolub
( Matrix<F>& A, Matrix<Int>& pPerm, Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    return BusingerGolub( A, t, d, pPerm, numSteps, alwaysRecompute );
}

template<typename F> 
inline Int
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, Matrix<Int>& pPerm, 
  bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    const Int numSteps = Min(A.Height(),A.Width());
    return BusingerGolub( A, t, d, pPerm, numSteps, alwaysRecompute );
}

// If we don't need 't' or 'd' from the above routine
template<typename F> 
inline Int
BusingerGolub( Matrix<F>& A, Matrix<Int>& pPerm, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    return BusingerGolub( A, t, d, pPerm, alwaysRecompute );
}

template<typename F>
inline ValueInt<Base<F>>
FindColPivot
( const DistMatrix<F>& A, const std::vector<Base<F>>& norms, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("qr::FindColPivot"))
    typedef Base<F> Real;
    const Int localColsBefore = A.LocalColOffset(col);
    const ValueInt<Real> localPivot = FindPivot( norms, localColsBefore );
    ValueInt<Real> pivot;
    pivot.value = localPivot.value;
    pivot.index = A.GlobalCol(localPivot.index);
    return mpi::AllReduce( pivot, mpi::MaxLocOp<Real>(), A.Grid().RowComm() );
}

template<typename F>
inline Base<F>
ColNorms( const DistMatrix<F>& A, std::vector<Base<F>>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ColNorms"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    mpi::Comm colComm = A.Grid().ColComm();
    mpi::Comm rowComm = A.Grid().RowComm();

    // Carefully perform the local portion of the computation
    std::vector<Real> localScales(localWidth,0), 
                      localScaledSquares(localWidth,1);
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));    
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= localScales[jLoc] )
                {
                    const Real relScale = alphaAbs/localScales[jLoc];
                    localScaledSquares[jLoc] += relScale*relScale;
                }
                else
                {
                    const Real relScale = localScales[jLoc]/alphaAbs;
                    localScaledSquares[jLoc] = 
                        localScaledSquares[jLoc]*relScale*relScale + 1;
                    localScales[jLoc] = alphaAbs;
                }
            }
        }
    }

    // Find the maximum relative scales 
    std::vector<Real> scales(localWidth);
    mpi::AllReduce
    ( localScales.data(), scales.data(), localWidth, mpi::MAX, colComm );

    // Equilibrate the local scaled sums to the maximum scale
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        if( scales[jLoc] != 0 )
        {
            const Real relScale = localScales[jLoc]/scales[jLoc];
            localScaledSquares[jLoc] *= relScale*relScale;
        }
    }

    // Now sum the local contributions (can ignore results where scale is 0)
    std::vector<Real> scaledSquares(localWidth); 
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), localWidth, colComm );

    // Finish the computation
    Real maxLocalNorm = 0;
    norms.resize( localWidth );
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        if( scales[jLoc] != 0 )
            norms[jLoc] = scales[jLoc]*Sqrt(scaledSquares[jLoc]);
        else
            norms[jLoc] = 0;
        maxLocalNorm = std::max( maxLocalNorm, norms[jLoc] );
    }
    return mpi::AllReduce( maxLocalNorm, mpi::MAX, rowComm );
}

template<typename F>
inline void
ReplaceColNorms
( const DistMatrix<F>& A, std::vector<Int>& inaccurateNorms, 
  std::vector<Base<F>>& norms, std::vector<Base<F>>& origNorms )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ReplaceColNorms"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const Int numInaccurate = inaccurateNorms.size();
    mpi::Comm colComm = A.Grid().ColComm();

    // Carefully perform the local portion of the computation
    std::vector<Real> localScales(numInaccurate,0), 
                      localScaledSquares(numInaccurate,1);
    for( Int s=0; s<numInaccurate; ++s )
    {
        const Int jLoc = inaccurateNorms[s];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));    
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= localScales[s] )
                {
                    const Real relScale = alphaAbs/localScales[s];
                    localScaledSquares[s] += relScale*relScale;
                }
                else
                {
                    const Real relScale = localScales[s]/alphaAbs;
                    localScaledSquares[s] = 
                        localScaledSquares[s]*relScale*relScale + 1;
                    localScales[s] = alphaAbs;
                }
            }
        }
    }

    // Find the maximum relative scales 
    std::vector<Real> scales(numInaccurate);
    mpi::AllReduce
    ( localScales.data(), scales.data(), numInaccurate, mpi::MAX, colComm );

    // Equilibrate the local scaled sums to the maximum scale
    for( Int s=0; s<numInaccurate; ++s )
    {
        if( scales[s] != 0 )
        {
            const Real relScale = localScales[s]/scales[s];
            localScaledSquares[s] *= relScale*relScale;
        }
    }

    // Now sum the local contributions (can ignore results where scale is 0)
    std::vector<Real> scaledSquares(numInaccurate); 
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), numInaccurate, colComm );

    // Finish the computation
    for( Int s=0; s<numInaccurate; ++s )
    {
        const Int jLoc = inaccurateNorms[s];
        if( scales[s] != 0 )
            norms[jLoc] = scales[s]*Sqrt(scaledSquares[s]);
        else
            norms[jLoc] = 0;
        origNorms[jLoc] = norms[jLoc];
    }
}

template<typename F,Dist UPerm>
inline Int
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d, 
  DistMatrix<Int,UPerm,STAR>& pPerm, Int maxSteps, Base<F> tol, 
  bool alwaysRecompute=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("qr::BusingerGolub");
        if( A.Grid() != pPerm.Grid() || A.Grid() != t.Grid() || 
            t.Grid() != d.Grid() )
            LogicError("A, t, d, and pPerm must have the same grid");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    maxSteps = Min(maxSteps,Min(m,n));
    t.SetRoot( A.DiagonalRoot() );
    d.SetRoot( A.DiagonalRoot() );
    t.AlignCols( A.DiagonalAlign() );
    d.AlignCols( A.DiagonalAlign() );
    t.Resize( maxSteps, 1 );
    d.Resize( maxSteps, 1 );

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms( A.LocalWidth() );
    const Real maxOrigNorm = ColNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());
    std::vector<Int> inaccurateNorms;

    // Initialize the inverse permutation to the identity
    DistMatrix<Int,UPerm,STAR> pInvPerm( pPerm.Grid() );
    pInvPerm.AlignWith( pPerm );
    pInvPerm.Resize( n, 1 );
    for( Int jLoc=0; jLoc<pInvPerm.LocalHeight(); ++jLoc ) 
        pInvPerm.SetLocal( jLoc, 0, pInvPerm.GlobalRow(jLoc) );

    const Grid& g = A.Grid();
    DistMatrix<F> z21(g);
    DistMatrix<F,MC,STAR> aB1_MC_STAR(g);
    DistMatrix<F,MR,STAR> z21_MR_STAR(g);
    DistMatrix<F,STAR,MR> a12_STAR_MR(g);

    Int k=0;
    for( ; k<maxSteps; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find the next column pivot
        const ValueInt<Real> pivot = FindColPivot( A, norms, k );
        if( pivot.value <= tol*maxOrigNorm )
            break;
        RowSwap( pInvPerm, k, pivot.index );

        // Perform the swap
        const Int jPiv = pivot.index;
        const Int curOwner = A.ColOwner(k);
        const Int pivOwner = A.ColOwner(jPiv);
        const Int myCur = A.IsLocalCol(k);
        const Int myPiv = A.IsLocalCol(jPiv);
        if( jPiv != k )
        {
            if( myCur && myPiv )
            {
                const Int kLoc    = A.LocalCol(k);
                const Int jPivLoc = A.LocalCol(jPiv);
                blas::Swap
                ( mLocal, A.Buffer(0,kLoc), 1, A.Buffer(0,jPivLoc), 1 );
                norms[jPivLoc] = norms[kLoc];
                origNorms[jPivLoc] = origNorms[kLoc];
            }
            else if( myCur )
            {
                const Int kLoc = A.LocalCol(k);
                mpi::SendRecv
                ( A.Buffer(0,kLoc), mLocal, pivOwner, pivOwner, g.RowComm() );
                mpi::Send( norms[kLoc], pivOwner, g.RowComm() );
            }
            else if( myPiv )
            {
                const Int jPivLoc = A.LocalCol(jPiv);
                mpi::SendRecv
                ( A.Buffer(0,jPivLoc), mLocal, 
                  curOwner, curOwner, g.RowComm() );
                norms[jPivLoc] = mpi::Recv<Real>( curOwner, g.RowComm() );
            }
        }

        // Find tau and u such that
        //  / I - tau | 1 | | 1, u^H | \ | alpha11 | = | beta |
        //  \         | u |            / |     a21 | = |    0 |
        const F tau = LeftReflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        F alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        // AB2 := Hous(aB1,tau) AB2
        //      = (I - tau aB1 aB1^H) AB2
        //      = AB2 - tau aB1 (AB2^H aB1)^H
        aB1_MC_STAR.AlignWith( AB2 );
        aB1_MC_STAR = aB1;
        z21_MR_STAR.AlignWith( AB2 );
        Zeros( z21_MR_STAR, AB2.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC_STAR, F(0), z21_MR_STAR );
        z21_MR_STAR.SumOver( AB2.ColComm() );
        Ger
        ( -tau, aB1_MC_STAR.LockedMatrix(), z21_MR_STAR.LockedMatrix(),
          AB2.Matrix() );

        // Reset alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176.
        // However, we do so in two steps in order to lower the communication
        // latency:
        //   1) Each process first computes which of its column norms are
        //      too inaccurate and need to be recomputed.
        //   2) Each process communicates within its process column in order
        //      to replace the inaccurate column norms.
        // Step 1: Perform all of the easy updates and mark inaccurate norms
        a12_STAR_MR = a12;
        inaccurateNorms.resize(0);
        const Int a12LocalWidth = a12_STAR_MR.LocalWidth();
        for( Int jLoc12=0; jLoc12<a12LocalWidth; ++jLoc12 )
        {
            const Int j = (k+1) + a12.GlobalCol(jLoc12);
            const Int jLoc = A.LocalCol(j);
            if( norms[jLoc] != Real(0) )
            {
                const Real beta = Abs(a12_STAR_MR.GetLocal(0,jLoc12));
                Real gamma = beta / norms[jLoc];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );

                const Real ratio = norms[jLoc] / origNorms[jLoc];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || alwaysRecompute )
                    inaccurateNorms.push_back( jLoc );
                else
                    norms[jLoc] *= Sqrt(gamma);
            }
        }
        // Step 2: Compute the replacement norms and also reset origNorms
        ReplaceColNorms( A, inaccurateNorms, norms, origNorms );
    }
    InvertPermutation( pInvPerm, pPerm );

    // Form d and rescale R
    auto R = View( A, 0, 0, k, n );
    d = R.GetRealPartOfDiagonal();
    const Int diagLengthLoc = d.LocalHeight();
    for( Int jLoc=0; jLoc<diagLengthLoc; ++jLoc )
    {
        const Real delta = d.GetLocal(jLoc,0);
        if( delta >= Real(0) )
            d.SetLocal(jLoc,0,Real(1));
        else
            d.SetLocal(jLoc,0,Real(-1));
    }
    DiagonalScaleTrapezoid( LEFT, UPPER, NORMAL, d, R );

    return k;
}

// If we don't need 't' or 'd' from the above routine
template<typename F,Dist UPerm>
inline Int
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm,
  Int maxSteps, Base<F> tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    DistMatrix<F,MD,STAR> t( A.Grid() );
    DistMatrix<Base<F>,MD,STAR> d( A.Grid() );
    return BusingerGolub( A, t, d, pPerm, maxSteps, tol, alwaysRecompute );
}

template<typename F,Dist UPerm>
inline Int
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d, 
  DistMatrix<Int,UPerm,STAR>& pPerm, Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    // Use a tolerance of -1 so that we do not stop early
    return BusingerGolub
           ( A, t, d, pPerm, numSteps, Base<F>(-1), alwaysRecompute );
}

// If we don't need 't' or 'd' from the above routine
template<typename F,Dist UPerm>
inline Int
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm,
  Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    DistMatrix<F,MD,STAR> t( A.Grid() );
    DistMatrix<Base<F>,MD,STAR> d( A.Grid() );
    return BusingerGolub( A, t, d, pPerm, numSteps, alwaysRecompute );
}

template<typename F,Dist UPerm>
inline Int
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d,
  DistMatrix<Int,UPerm,STAR>& pPerm, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    const Int numSteps = Min(A.Height(),A.Width());
    return BusingerGolub( A, t, d, pPerm, numSteps, alwaysRecompute );
}

// If we don't need 't' or 'd' from the above routine
template<typename F,Dist UPerm>
inline Int
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm, 
  bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    DistMatrix<F,MD,STAR> t( A.Grid() );
    DistMatrix<Base<F>,MD,STAR> d( A.Grid() );
    return BusingerGolub( A, t, d, pPerm, alwaysRecompute );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_QR_BUSINGERGOLUB_HPP
