/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_APPLYPACKEDREFLECTORS_RLVB_HPP
#define ELEM_APPLYPACKEDREFLECTORS_RLVB_HPP

#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SETDIAGONAL_INC
#include ELEM_GEMM_INC
#include ELEM_SYRK_INC
#include ELEM_HERK_INC
#include ELEM_TRSM_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace apply_packed_reflectors {

//
// Since applying Householder transforms from vectors stored right-to-left
// implies that we will be forming a generalization of
//
//   (I - tau_1 u_1 u_1^H) (I - tau_0 u_0 u_0^H) = 
//   I - tau_0 u_0 u_0^H - tau_1 u_1 u_1^H + (tau_0 tau_1 u_1^H u_0) u_1 u_0^H =
//   I - [ u_0, u_1 ] [  tau_0,                 0     ] [ u_0^H ]
//                    [ -tau_0 tau_1 u_1^H u_0, tau_1 ] [ u_1^H ],
//
// which has a lower-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   tril(T) = tril( U^H U ),  diag(T) = 1/t or 1/conj(t),
//
// where U is the matrix of Householder vectors and t is the vector of scalars.
//

template<typename F> 
inline void
RLVB
( Conjugation conjugation, Int offset, 
  const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("apply_packed_reflectors::RLVB");
        if( H.Height() != A.Width() )
            LogicError("H's height must match A's width");
    )
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    DEBUG_ONLY(
        if( t.Height() != diagLength )
            LogicError("t must be the same length as H's offset diag");
    )
    Matrix<F> HPanCopy, SInv, Z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( diagLength, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan = LockedViewRange( H, ki, kj, nA, kj+nb );
        auto ARight = ViewRange( A, 0, ki, mA, nA );
        auto t1 = LockedView( t, k, 0, nb, 1 );

        HPanCopy = HPan;
        MakeTriangular( LOWER, HPanCopy );
        SetDiagonal( HPanCopy, F(1) );

        Herk( LOWER, ADJOINT, F(1), HPanCopy, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, NORMAL, F(1), ARight, HPanCopy, Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( NORMAL, ADJOINT, F(-1), Z, HPanCopy, F(1), ARight );
    }
}

template<typename F> 
inline void
RLVB
( Conjugation conjugation, Int offset, 
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("apply_packed_reflectors::RLVB");
        if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
            LogicError("{H,t,A} must be distributed over the same grid");
    )
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    DEBUG_ONLY(
        if( t.Height() != diagLength )
            LogicError("t must be the same length as H's offset diag");
        if( !H.DiagonalAlignedWith( t, offset ) )
            LogicError("t must be aligned with H's 'offset' diagonal");
    )
    const Grid& g = H.Grid();
    DistMatrix<F> HPanCopy(g);
    DistMatrix<F,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<F,MR,  STAR> HPan_MR_STAR(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > ZAdj_STAR_MC(g);
    DistMatrix<F,STAR,VC  > ZAdj_STAR_VC(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( diagLength, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan = LockedViewRange( H, ki, kj, nA, kj+nb );
        auto ARight = ViewRange( A, 0, ki, mA, nA );
        auto t1 = LockedView( t, k, 0, nb, 1 );

        HPanCopy = HPan;
        MakeTriangular( LOWER, HPanCopy );
        SetDiagonal( HPanCopy, F(1) );

        HPan_VC_STAR = HPanCopy;
        Zeros( SInv_STAR_STAR, nb, nb );
        Herk
        ( LOWER, ADJOINT, 
          F(1), HPan_VC_STAR.LockedMatrix(),
          F(0), SInv_STAR_STAR.Matrix() );     
        SInv_STAR_STAR.SumOver( HPan_VC_STAR.ColComm() );
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MR_STAR.AlignWith( ARight );
        HPan_MR_STAR = HPan_VC_STAR;
        ZAdj_STAR_MC.AlignWith( ARight );
        LocalGemm( ADJOINT, ADJOINT, F(1), HPan_MR_STAR, ARight, ZAdj_STAR_MC );
        ZAdj_STAR_VC.AlignWith( ARight );
        ZAdj_STAR_VC.PartialRowSumScatterFrom( ZAdj_STAR_MC );
 
        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, ADJOINT, F(-1), ZAdj_STAR_MC, HPan_MR_STAR, F(1), ARight );
    }
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef ELEM_APPLYPACKEDREFLECTORS_RLVB_HPP
