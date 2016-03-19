/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTQUASITRSM_LLT_HPP
#define ELEM_MULTISHIFTQUASITRSM_LLT_HPP

#include ELEM_GEMM_INC
#include ELEM_TRSM_INC

namespace elem {
namespace msquasitrsm {

template<typename F>
inline void
LLTUnb
( bool conjugate, const Matrix<F>& L, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("msquasitrsm::LLTUnb"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    const F* LBuf = L.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldl = L.LDim();
    const Int ldx = X.LDim();

    if( conjugate )
        Conjugate( X );

    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && LBuf[(k-1)+k*ldl] != F(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear systems via 2x2 LQ decompositions produced
            // by the Givens rotation
            //    | L(k,k)-shift L(k,k+1) | | c -conj(s) | = | gamma11 0 |
            //                              | s    c     |
            // and by also forming the bottom two entries of the 2x2 resulting
            // lower-triangular matrix, say gamma21 and gamma22
            //
            // Extract the constant part of the 2x2 diagonal block, D
            const F delta12 = LBuf[   k +(k+1)*ldl];
            const F delta21 = LBuf[(k+1)+   k *ldl];
            for( Int j=0; j<n; ++j )
            {
                const F delta11 = LBuf[   k +   k *ldl] - shifts.Get(j,0);
                const F delta22 = LBuf[(k+1)+(k+1)*ldl] - shifts.Get(j,0);
                // Decompose D = L Q
                Real c; F s;
                const F gamma11 = blas::Givens( delta11, delta12, &c, &s );
                const F gamma21 =        c*delta21 + s*delta22;
                const F gamma22 = -Conj(s)*delta21 + c*delta22;

                F* xBuf = &XBuf[j*ldx];
                // Solve against Q^T
                const F chi1 = xBuf[k  ];
                const F chi2 = xBuf[k+1];
                xBuf[k  ] =        c*chi1 + s*chi2;
                xBuf[k+1] = -Conj(s)*chi1 + c*chi2;

                // Solve against R^T
                xBuf[k+1] /= gamma22;
                xBuf[k  ] -= gamma21*xBuf[k+1];
                xBuf[k  ] /= gamma11;

                // Update x0 := x0 - L10^T x1
                blas::Axpy( k, -xBuf[k  ], &LBuf[k  ], ldl, xBuf, 1 );
                blas::Axpy( k, -xBuf[k+1], &LBuf[k+1], ldl, xBuf, 1 );
            }
        }
        else
        {
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                // Solve the 1x1 linear system
                xBuf[k] /= LBuf[k+k*ldl] - shifts.Get(j,0);

                // Update x0 := x0 - l10^T chi_1
                blas::Axpy( k, -xBuf[k], &LBuf[k], ldl, xBuf, 1 );
            }
        }
        --k;
    }
    if( conjugate )
        Conjugate( X );
}

template<typename F>
inline void
LLT
( Orientation orientation, 
  const Matrix<F>& L, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LLT");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();

    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( X );

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        LLTUnb( false, L11, shifts, X1 );
        Gemm( TRANSPOSE, NORMAL, F(-1), L10, X1, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }

    if( conjugate )
        Conjugate( X );
}

// width(X) >> p
template<typename F>
inline void
LLTLarge
( Orientation orientation, 
  const DistMatrix<F>& L, const DistMatrix<F,VR,STAR>& shifts, 
        DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LLTLarge");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-[T/H][* ,* ] X1[* ,VR]
        LocalMultiShiftQuasiTrsm
        ( LEFT, LOWER, orientation, F(1), L11_STAR_STAR, shifts, X1_STAR_VR );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR] <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;        // L10[* ,MC] <- L10[MC,MR]

        // X0[MC,MR] -= (L10[* ,MC])^(T/H) X1[* ,MR]
        //            = L10^[T/H][MC,* ] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), L10_STAR_MC, X1_STAR_MR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

// width(X) ~= p
template<typename F,Dist shiftColDist,Dist shiftRowDist>
inline void
LLTMedium
( Orientation orientation, const DistMatrix<F>& L, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts, DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LLTMedium");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    DistMatrix<F,MR,  STAR> shifts_MR_STAR( shifts ),
                            shifts_MR_STAR_Align(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        // X1[* ,MR] <- X1[MC,MR]
        X1Trans_MR_STAR.AlignWith( X0 );
        X1.TransposeColAllGather( X1Trans_MR_STAR, (orientation==ADJOINT) );

        // X1[* ,MR] := L11^-[T/H][* ,* ] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] L11^-1[* ,* ]
        shifts_MR_STAR_Align.AlignWith( X1Trans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR;
        LocalMultiShiftQuasiTrsm
        ( RIGHT, LOWER, NORMAL,
          F(1), L11_STAR_STAR, shifts_MR_STAR_Align, X1Trans_MR_STAR );

        X1.TransposeColFilterFrom( X1Trans_MR_STAR, (orientation==ADJOINT) );
        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10; // L10[* ,MC] <- L10[MC,MR]

        // X0[MC,MR] -= (L10[* ,MC])^[T/H] X1[* ,MR]
        //            = L10^[T/H][MC,* ] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), L10_STAR_MC, X1Trans_MR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

using trsm::AddInLocalData;

// width(X) << p
template<typename F,Dist colDist,Dist shiftColDist,Dist shiftRowDist>
inline void
LLTSmall
( Orientation orientation, 
  const DistMatrix<F,colDist,STAR>& L, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts, 
        DistMatrix<F,colDist,STAR>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LLTSmall");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( L.Height() != L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X must be aligned");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), Z1_STAR_STAR(g),
                            shifts_STAR_STAR(shifts);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        // X1 -= L21' X2
        LocalGemm( orientation, NORMAL, F(-1), L21, X2, Z1_STAR_STAR );
        AddInLocalData( X1, Z1_STAR_STAR );
        Z1_STAR_STAR.SumOver( X1.DistComm() );

        // X1 := L11^-1 X1
        L11_STAR_STAR = L11;
        LocalMultiShiftQuasiTrsm
        ( LEFT, LOWER, orientation, F(1), 
          L11_STAR_STAR, shifts_STAR_STAR, Z1_STAR_STAR );
        X1 = Z1_STAR_STAR;

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F,Dist rowDist,Dist shiftColDist,Dist shiftRowDist>
inline void
LLTSmall
( Orientation orientation, 
  const DistMatrix<F,STAR,rowDist>& L, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts, 
        DistMatrix<F,rowDist,STAR>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LLTSmall");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( L.Height() != L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
        if( L.RowAlign() != X.ColAlign() )
            LogicError("L and X must be aligned");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X1_STAR_STAR(g),
                            shifts_STAR_STAR(shifts);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[* ,VR]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VR,* ]

        // X1[* ,* ] := L11^-[T/H][* ,* ] X1[* ,* ]
        LocalMultiShiftQuasiTrsm
        ( LEFT, LOWER, orientation,
          F(1), L11_STAR_STAR, shifts_STAR_STAR, X1_STAR_STAR );

        X1 = X1_STAR_STAR;

        // X0[VR,* ] -= L10[* ,VR]^(T/H) X1[* ,* ]
        LocalGemm( orientation, NORMAL, F(-1), L10, X1_STAR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

} // namespace msquasitrsm
} // namespace elem

#endif // ifndef ELEM_MULTISHIFTQUASITRSM_LLT_HPP
