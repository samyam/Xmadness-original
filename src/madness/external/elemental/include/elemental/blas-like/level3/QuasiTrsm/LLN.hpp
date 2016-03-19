/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QUASITRSM_LLN_HPP
#define ELEM_QUASITRSM_LLN_HPP

#include ELEM_GEMM_INC

namespace elem {
namespace quasitrsm {

// Left Lower NORMAL (Non)Unit QuasiTrsm 
//   X := tril(L)^-1  X, or
//   X := trilu(L)^-1 X

template<typename F>
inline void
LLNUnb( const Matrix<F>& L, Matrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("quasitrsm::LLNUnb"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    const F* LBuf = L.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldl = L.LDim();
    const Int ldx = X.LDim();

    Int k=0;
    while( k < m )
    {
        const bool in2x2 = ( k+1<m && LBuf[k+(k+1)*ldl] != F(0) );
        if( in2x2 ) 
        {
            // Solve the 2x2 linear systems via a 2x2 LQ decomposition produced
            // by the Givens rotation
            //    | L(k,k) L(k,k+1) | | c -conj(s) | = | gamma11 0 |
            //                        | s    c     |
            // and by also forming the bottom two entries of the 2x2 resulting
            // lower-triangular matrix, say gamma21 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = LBuf[ k   + k   *ldl];
            const F delta12 = LBuf[ k   +(k+1)*ldl];
            const F delta21 = LBuf[(k+1)+ k   *ldl];
            const F delta22 = LBuf[(k+1)+(k+1)*ldl];
            // Decompose D = L Q
            Real c; F s;
            const F gamma11 = lapack::Givens( delta11, delta12, &c, &s );
            const F gamma21 =        c*delta21 + s*delta22;
            const F gamma22 = -Conj(s)*delta21 + c*delta22;
            if( checkIfSingular )
            {
                // TODO: Check if sufficiently small instead
                if( gamma11 == F(0) || gamma22 == F(0) )
                    LogicError("Singular diagonal block detected");
            }
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];

                // Solve against L
                xBuf[k  ] /= gamma11;
                xBuf[k+1] -= gamma21*xBuf[k];
                xBuf[k+1] /= gamma22;

                // Solve against Q
                const F chi1 = xBuf[k  ];
                const F chi2 = xBuf[k+1];
                xBuf[k  ] = c*chi1 - Conj(s)*chi2;
                xBuf[k+1] = s*chi1 +       c*chi2;

                // Update x2 := x2 - L21 x1
                blas::Axpy
                ( m-(k+2), -xBuf[k  ], 
                  &LBuf[(k+2)+ k   *ldl], 1, &xBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xBuf[k+1], 
                  &LBuf[(k+2)+(k+1)*ldl], 1, &xBuf[k+2], 1 );
            }

            k += 2;
        }
        else
        {
            if( checkIfSingular )
            {
                // TODO: Check if sufficiently small instead
                if( LBuf[k+k*ldl] == F(0) )
                    LogicError("Singular diagonal entry detected");
            }
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                xBuf[k] /= LBuf[k+k*ldl];
                blas::Axpy
                ( m-(k+1), -xBuf[k], &LBuf[(k+1)+k*ldl], 1, &xBuf[k+1], 1 );
            }
            k += 1;
        }
    }
}

template<typename F>
inline void
LLN( const Matrix<F>& L, Matrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("quasitrsm::LLN"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        LLNUnb( L11, X1, checkIfSingular );
        Gemm( NORMAL, NORMAL, F(-1), L21, X1, F(1), X2 );
    }
}

// For large numbers of RHS's, e.g., width(X) >> p
template<typename F>
inline void
LLNLarge( const DistMatrix<F>& L, DistMatrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("quasitrsm::LLNLarge"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-1[* ,* ] X1[* ,VR]
        LocalQuasiTrsm
        ( LEFT, LOWER, NORMAL, F(1), L11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;        // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), L21_MC_STAR, X1_STAR_MR, F(1), X2 );
    }
}

// For medium numbers of RHS's, e.g., width(X) ~= p
template<typename F>
inline void
LLNMedium( const DistMatrix<F>& L, DistMatrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("quasitrsm::LLNMedium"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1Trans_MR_STAR.AlignWith( X2 );
        X1.TransposeColAllGather( X1Trans_MR_STAR ); // X1[* ,MR] <- X1[MC,MR]

        // X1^T[MR,* ] := X1^T[MR,* ] L11^-T[* ,* ]
        //              = (L11^-1[* ,* ] X1[* ,MR])^T
        LocalQuasiTrsm
        ( RIGHT, LOWER, TRANSPOSE,
          F(1), L11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        X1.TransposeColFilterFrom( X1Trans_MR_STAR ); // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;                   // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, F(-1), L21_MC_STAR, X1Trans_MR_STAR, F(1), X2 );
    }
}

// For small numbers of RHS's, e.g., width(X) < p
template<typename F,Dist colDist>
inline void
LLNSmall
( const DistMatrix<F,colDist,STAR>& L, DistMatrix<F,colDist,STAR>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("quasitrsm::LLNSmall");
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X1_STAR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        LocalQuasiTrsm
        ( LEFT, LOWER, NORMAL, 
          F(1), L11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
    }
}

} // namespace quasitrsm
} // namespace elem

#endif // ifndef ELEM_QUASITRSM_LLN_HPP
