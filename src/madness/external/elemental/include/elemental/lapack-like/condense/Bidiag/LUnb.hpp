/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BIDIAG_LUNB_HPP
#define ELEM_BIDIAG_LUNB_HPP

#include ELEM_CONJUGATE_INC
#include ELEM_GEMV_INC
#include ELEM_GER_INC

#include ELEM_REFLECTOR_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace bidiag {

template<typename F>
inline void LUnb( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::LUnb");
        if( A.Height() > A.Width() )
            LogicError("A must be at least as wide as it is tall");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int tPHeight = m;
    const Int tQHeight = Max(m-1,0);
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    // Views
    Matrix<F> alpha21T, a21B;

    // Temporaries
    Matrix<F> x12Adj, w21;

    for( Int k=0; k<m; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );
        auto a1R     = ViewRange( A, k,   k,   k+1, n   );
        auto A2R     = ViewRange( A, k+1, k,   m,   n   );

        // Find tauP and v such that
        //  | alpha11 a12 | / I - tauP | 1   | | 1, conj(v) | \ = | epsilonP 0 |
        //                  \          | v^T |                /
        const F tauP = RightReflector( alpha11, a12 );
        tP.Set(k,0,tauP);

        // Temporarily set a1R = | 1 v |
        const F epsilonP = alpha11.Get(0,0);
        alpha11.Set(0,0,F(1));

        // A2R := A2R Hous(a1R^T,tauP)
        //      = A2R (I - tauP a1R^T conj(a1R))
        //      = A2R - tauP (A2R a1R^T) conj(a1R)
        // ---------------------------------------
        // w21 := A2R a1R^T = A2R |   1 |
        //                        | v^T |
        Zeros( w21, a21.Height(), 1 );
        Gemv( NORMAL, F(1), A2R, a1R, F(0), w21 );
        // A2R := A2R - tauP w21 conj(a1R)
        Ger( -tauP, w21, a1R, A2R );

        // Put epsilonP back 
        alpha11.Set(0,0,epsilonP);

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown( a21, alpha21T, a21B, 1 );

            // Find tauQ and u such that
            //  / I - tauQ | 1 | | 1, u^H | \ | alpha21T | = | epsilonQ |
            //  \          | u |            / | a21B     | = |    0     |
            const F tauQ = LeftReflector( alpha21T, a21B );
            tQ.Set(k,0,tauQ);

            // Temporarily set a21 = | 1 |
            //                       | u |
            const F epsilonQ = alpha21T.Get(0,0);
            alpha21T.Set(0,0,F(1));

            // A22 := Hous(a21,tauQ) A22
            //      = (I - tauQ a21 a21^H) A22
            //      = A22 - tauQ a21 (A22^H a21)^H
            // -----------------------------------
            // x12^H := (a21^H A22)^H = A22^H a21
            Zeros( x12Adj, a12.Width(), 1 );
            Gemv( ADJOINT, F(1), A22, a21, F(0), x12Adj );
            // A22 := A22 - tauQ a21 x12 
            //      = (I - tauQ a21 a21^H) A22
            Ger( -tauQ, a21, x12Adj, A22 );

            // Put epsilonQ back
            alpha21T.Set(0,0,epsilonQ);
        }
    }
}

template<typename F> 
inline void LUnb
( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::LUnb");
        if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
            LogicError("Process grids do not match");
        if( A.Height() > A.Width() )
            LogicError("A must be at least as wide as it is tall");
    )
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int tPHeight = m;
    const Int tQHeight = Max(m-1,0);
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    DistMatrix<F,MC,  STAR> a21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > a1R_STAR_MR(g);
    DistMatrix<F,MR,  STAR> x12Adj_MR_STAR(g);
    DistMatrix<F,MC,  STAR> w21_MC_STAR(g);

    for( Int k=0; k<m; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );
        auto a1R     = ViewRange( A, k,   k,   k+1, n   );
        auto A2R     = ViewRange( A, k+1, k,   m,   n   );

        // Find tauP and v such that
        //  | alpha11 a12 | / I - tauP | 1   | | 1, conj(v) | \ = | epsilonP 0 |
        //                  \          | v^T |                /
        const F tauP = RightReflector( alpha11, a12 );
        tP.Set(k,0,tauP);

        // Temporarily set a1R = | 1 v |
        F epsilonP=0;
        if( alpha11.IsLocal(0,0) )
            epsilonP = alpha11.GetLocal(0,0);
        alpha11.Set(0,0,F(1));

        // A2R := A2R Hous(a1R^T,tauP)
        //      = A2R (I - tauP a1R^T conj(a1R))
        //      = A2R - tauP (A2R a1R^T) conj(a1R)
        // -------------------------------------
        // w21 := A2R a1R^T = A2R | 1   |
        //                        | v^T |
        alpha11.Set(0,0,F(1));
        a1R_STAR_MR.AlignWith( A2R );
        a1R_STAR_MR = a1R;
        w21_MC_STAR.AlignWith( A2R );
        Zeros( w21_MC_STAR, a21.Height(), 1 );
        LocalGemv( NORMAL, F(1), A2R, a1R_STAR_MR, F(0), w21_MC_STAR );
        w21_MC_STAR.SumOver( A2R.RowComm() );
        // A2R := A2R - tauP w21 conj(a1R)
        LocalGer( -tauP, w21_MC_STAR, a1R_STAR_MR, A2R );

        // Put epsilonP back 
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,epsilonP);

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            DistMatrix<F> alpha21T(g), a21B(g);
            PartitionDown( a21, alpha21T, a21B, 1 );

            // Find tauQ and u such that
            //  / I - tauQ | 1 | | 1, u^H | \ | alpha21T | = | epsilonQ |
            //  \          | u |            / | a21B     | = |    0     |
            const F tauQ = LeftReflector( alpha21T, a21B );
            tQ.Set(k,0,tauQ);

            // Temporarily set a21 = | 1 |
            //                       | u |
            F epsilonQ=0;
            if( alpha21T.IsLocal(0,0) )
                epsilonQ = alpha21T.GetLocal(0,0);
            alpha21T.Set(0,0,F(1));

            // A22 := Hous(a21,tauQ) A22
            //      = (I - tauQ a21 a21^H) A22
            //      = A22 - tauQ a21 (A22^H a21)^H
            // ----------------------------------
            // x12^H := (a21^H A22)^H = A22^H a21
            a21_MC_STAR.AlignWith( A22 );
            a21_MC_STAR = a21;
            x12Adj_MR_STAR.AlignWith( A22 );
            Zeros( x12Adj_MR_STAR, a12.Width(), 1 );
            LocalGemv( ADJOINT, F(1), A22, a21_MC_STAR, F(0), x12Adj_MR_STAR );
            x12Adj_MR_STAR.SumOver( A22.ColComm() );
            // A22 := A22 - tauQ a21 x12
            //      = (I - tauQ a21 a21^H) A22
            LocalGer( -tauQ, a21_MC_STAR, x12Adj_MR_STAR, A22 );

            // Put epsilonQ back
            if( alpha21T.IsLocal(0,0) )
                alpha21T.SetLocal(0,0,epsilonQ);
        }
    }
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_BIDIAG_LUNB_HPP
