/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BIDIAG_UUNB_HPP
#define ELEM_BIDIAG_UUNB_HPP

#include ELEM_CONJUGATE_INC
#include ELEM_GEMV_INC
#include ELEM_GER_INC

#include ELEM_REFLECTOR_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace bidiag {

template<typename F>
inline void UUnb( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::UUnb"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
    DEBUG_ONLY(
        if( m < n )
            LogicError("A must be at least as tall as it is wide");
    )
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    Matrix<F> x12Adj, w21;

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha11 | = | epsilonQ |
        //  \          | u |            / |     a21 | = |    0     |
        const F tauQ = LeftReflector( alpha11, a21 );
        tQ.Set(k,0,tauQ );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        const F epsilonQ = alpha11.Get(0,0);
        alpha11.Set(0,0,F(1));

        // AB2 := Hous(aB1,tauQ) AB2
        //      = (I - tauQ aB1 aB1^H) AB2
        //      = AB2 - tauQ aB1 (AB2^H aB1)^H
        // -----------------------------------
        // x12^H := (aB1^H AB2)^H = AB2^H aB1
        Zeros( x12Adj, a12.Width(), 1 );
        Gemv( ADJOINT, F(1), AB2, aB1, F(0), x12Adj );
        // AB2 := AB2 - tauQ aB1 x12
        //      = (I - tauQ aB1 aB1^H) AB2
        Ger( -tauQ, aB1, x12Adj, AB2 );

        // Put epsilonQ back 
        alpha11.Set(0,0,epsilonQ);

        if( k+1 < n )
        {
            // Expose the subvector we seek to zero, a12R
            Matrix<F> alpha12L, a12R;
            PartitionRight( a12, alpha12L, a12R, 1 );

            // Find tauP and v such that
            //  |alpha12L a12R| /I - tauP |1  | |1, conj(v)|\ = |epsilonP 0|
            //                  \         |v^T|             /
            const F tauP = RightReflector( alpha12L, a12R );
            tP.Set(k,0,tauP);

            // Temporarily set a12^T = | 1 | 
            //                         | v |
            const F epsilonP = alpha12L.Get(0,0);
            alpha12L.Set(0,0,F(1));

            // A22 := A22 Hous(a12^T,tauP)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 - tauP (A22 a12^T) conj(a12)
            // ---------------------------------------
            // w21 := A22 a12^T = A22 | 1   |
            //                        | v^T |
            Zeros( w21, a21.Height(), 1 );
            Gemv( NORMAL, F(1), A22, a12, F(0), w21 );
            // A22 := A22 - tauP w21 conj(a12)
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back 
            alpha12L.Set(0,0,epsilonP);
        }
    }
}

template<typename F> 
inline void UUnb
( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::UUnb");
        if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
            LogicError("Process grids do not match");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( m < n )
            LogicError("A must be at least as tall as it is wide");
    )
    const Grid& g = A.Grid();
    const Int tPHeight = Max(n-1,0);
    const Int tQHeight = n;
    tP.Resize( tPHeight, 1 );
    tQ.Resize( tQHeight, 1 );

    DistMatrix<F,STAR,MR  > a12_STAR_MR(g);
    DistMatrix<F,MC,  STAR> aB1_MC_STAR(g);
    DistMatrix<F,MR,  STAR> x12Adj_MR_STAR(g);
    DistMatrix<F,MC,  STAR> w21_MC_STAR(g);

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha11 | = | epsilonQ |
        //  \          | u |            / |     a21 | = |    0     |
        const F tauQ = LeftReflector( alpha11, a21 );
        tQ.Set(k,0,tauQ );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        F epsilonQ=0;
        if( alpha11.IsLocal(0,0) )
            epsilonQ = alpha11.GetLocal(0,0);
        alpha11.Set(0,0,F(1));

        // AB2 := Hous(aB1,tauQ) AB2
        //      = (I - tauQ aB1 aB1^H) AB2
        //      = AB2 - tauQ aB1 (AB2^H aB1)^H
        // ----------------------------------
        // x12^H := (aB1^H AB2)^H = AB2^H aB1
        aB1_MC_STAR.AlignWith( aB1 );
        aB1_MC_STAR = aB1;
        x12Adj_MR_STAR.AlignWith( AB2 );
        Zeros( x12Adj_MR_STAR, a12.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC_STAR, F(0), x12Adj_MR_STAR );
        x12Adj_MR_STAR.SumOver( AB2.ColComm() );
        // AB2 := AB2 - tauQ aB1 x12
        LocalGer( -tauQ, aB1_MC_STAR, x12Adj_MR_STAR, AB2 );

        // Put epsilonQ back 
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,epsilonQ);

        if( k+1 < n )
        {
            // Expose the subvector we seek to zero, a12R
            DistMatrix<F> alpha12L(g), a12R(g);
            PartitionRight( a12, alpha12L, a12R, 1 );

            // Find tauP and v such that
            //  |alpha12L a12R| /I - tauP |1  | |1, conj(v)|\ = |epsilonP 0|
            //                  \         |v^T|             /
            const F tauP = RightReflector( alpha12L, a12R );
            tP.Set(k,0,tauP);

            // Temporarily set a12^T = | 1   |
            //                         | v^T |
            F epsilonP=0;
            if( alpha12L.IsLocal(0,0) )
                epsilonP = alpha12L.GetLocal(0,0);
            alpha12L.Set(0,0,F(1));

            // A22 := A22 Hous(a12^T,tauP)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 - tauP (A22 a12^T) conj(a12)
            // -------------------------------------
            // w21 := A22 a12^T = A22 | 1   |
            //                        | v^T |
            a12_STAR_MR.AlignWith( a12 );
            a12_STAR_MR = a12;
            w21_MC_STAR.AlignWith( A22 );
            Zeros( w21_MC_STAR, a21.Height(), 1 );
            LocalGemv( NORMAL, F(1), A22, a12_STAR_MR, F(0), w21_MC_STAR );
            w21_MC_STAR.SumOver( A22.RowComm() );
            // A22 := A22 - tauP w21 conj(a12)
            LocalGer( -tauP, w21_MC_STAR, a12_STAR_MR, A22 );

            // Put epsilonP back
            if( alpha12L.IsLocal(0,0) )
                alpha12L.SetLocal(0,0,epsilonP);
        }
    }
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_BIDIAG_UUNB_HPP
