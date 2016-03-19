/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QUASITRSV_UT_HPP
#define ELEM_QUASITRSV_UT_HPP

#include ELEM_AXPY_INC
#include ELEM_ZEROS_INC
#include ELEM_GEMV_INC

namespace elem {
namespace internal {

template<typename F>
inline void
QuasiTrsvUTUnb
( Orientation orientation, const Matrix<F>& U, Matrix<F>& x, 
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::QuasiTrsvUTUnb");
        if( U.Height() != U.Width() )
            LogicError("L must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal QuasiTrsvUT");
        if( orientation == NORMAL )
            LogicError("Invalid orientation");
    )
    typedef Base<F> Real;
    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( x );

    F* xBuf = x.Buffer();
    const F* UBuf = U.LockedBuffer();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int ldu = U.LDim();
    const Int m = U.Height();
    for( Int k=0; k<m; ++k )
    {
        const bool in2x2 = ( k+1<m && UBuf[(k+1)+k*ldu] != F(0) );
        if( in2x2 )
        {
            // Solve the 2x2 linear system via a 2x2 QR decomposition produced
            // by the Givens rotation
            //    | c        s | | U(k,  k) | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k) |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = UBuf[ k   + k   *ldu];
            const F delta12 = UBuf[ k   +(k+1)*ldu];
            const F delta21 = UBuf[(k+1)+ k   *ldu];
            const F delta22 = UBuf[(k+1)+(k+1)*ldu];
            // Decompose D = Q R
            Real c; F s;
            const F gamma11 = lapack::Givens( delta11, delta21, &c, &s );
            const F gamma12 =        c*delta12 + s*delta22;
            const F gamma22 = -Conj(s)*delta12 + c*delta22;
            if( checkIfSingular )
            {
                // TODO: Instead check if values are too small in magnitude
                if( gamma11 == F(0) || gamma22 == F(0) )
                    LogicError("Singular diagonal block detected");
            }
            // Solve against R^T
            xBuf[ k   *incx] /= gamma11;
            xBuf[(k+1)*incx] -= gamma12*xBuf[k*incx];
            xBuf[(k+1)*incx] /= gamma22;
            // Solve against Q^T
            const F chi1 = xBuf[ k   *incx];
            const F chi2 = xBuf[(k+1)*incx];
            xBuf[ k   *incx] = c*chi1 - Conj(s)*chi2;
            xBuf[(k+1)*incx] = s*chi1 +       c*chi2;

            // Update x2 := x2 - U12^T x1
            blas::Axpy
            ( m-(k+2), -xBuf[ k   *incx],
              &UBuf[ k   +(k+2)*ldu], ldu, &xBuf[(k+2)*incx], incx );
            blas::Axpy
            ( m-(k+2), -xBuf[(k+1)*incx],
              &UBuf[(k+1)+(k+2)*ldu], ldu, &xBuf[(k+2)*incx], incx );

            k += 2;
        }
        else
        {
            if( checkIfSingular )
                if( UBuf[k+k*ldu] == F(0) )
                    LogicError("Singular diagonal entry detected");
            // Solve the 1x1 linear system
            xBuf[k] /= UBuf[k+k*ldu];

            // Update x2 := x2 - u12^T chi_1
            blas::Axpy
            ( m-(k+1), -xBuf[k*incx],
              &UBuf[k+(k+1)*ldu], ldu, &xBuf[(k+1)*incx], incx );

            k += 1;
        }
        --k;
    }
    if( conjugate )
        Conjugate( x );
}

template<typename F>
inline void
QuasiTrsvUT
( Orientation orientation, const Matrix<F>& U, Matrix<F>& x,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::QuasiTrsvUT");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal QuasiTrsvUT");
        if( orientation == NORMAL )
            LogicError("Invalid orientation");
    )
    const bool vert = ( x.Width()==1 );
    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( x );

    Matrix<F> x1, x2;
    const Int m = U.Height();
    const Int bsize = Blocksize();
    Int k=0;
    while( k < m )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        if( vert )
        {
            x1 = ViewRange( x, k,    0, k+nb, 1 );
            x2 = ViewRange( x, k+nb, 0, m,    1 );
        }
        else
        {
            x1 = ViewRange( x, 0, k,    1, k+nb );
            x2 = ViewRange( x, 0, k+nb, 1, m    );
        }

        QuasiTrsvUTUnb( TRANSPOSE, U11, x1, checkIfSingular );
        Gemv( TRANSPOSE, F(-1), U12, x1, F(1), x2 );

        k += nb;
    }
    if( conjugate )
        Conjugate( x );
}

template<typename F>
inline void
QuasiTrsvUT
( Orientation orientation, const DistMatrix<F>& U, DistMatrix<F>& x,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::QuasiTrsvUT");
        if( U.Grid() != x.Grid() )
            LogicError("{U,x} must be distributed over the same grid");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal QuasiTrsvUT");
        if( orientation == NORMAL )
            LogicError("Invalid orientation");
    )
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();
    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( x );

    // Matrix views 
    DistMatrix<F> U11(g), U12(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MR,STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z1_MC_STAR(g), z2_MC_STAR(g);

        z_MC_STAR.AlignWith( U );
        Zeros( z_MC_STAR, m, 1 );

        Int k=0;
        while( k < m )
        {
            const Int nbProp = Min(bsize,m-k);
            const bool in2x2 =
                ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
            const Int nb = ( in2x2 ? nbProp+1 : nbProp );

            LockedViewRange( U11, U, k, k,    k+nb, k+nb );
            LockedViewRange( U12, U, k, k+nb, k+nb, m    );

            ViewRange( x1, x, k,    0, k+nb, 1 );

            ViewRange( z1_MC_STAR, z_MC_STAR, k,    0, k+nb, 1 );
            ViewRange( z2_MC_STAR, z_MC_STAR, k+nb, 0, m,    1 );

            if( k != 0 )
                x1.RowSumScatterUpdate( F(1), z1_MC_STAR );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            QuasiTrsvUT
            ( TRANSPOSE, U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(),
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_MR_STAR.AlignWith( U12 );
            x1_MR_STAR = x1_STAR_STAR;
            LocalGemv( TRANSPOSE, F(-1), U12, x1_MR_STAR, F(1), z2_MC_STAR );

            k += nb;
        }
    }
    else
    {
        DistMatrix<F,STAR,MR  > x1_STAR_MR(g);
        DistMatrix<F,MC,  MR  > z1(g);
        DistMatrix<F,MR,  MC  > z1_MR_MC(g);
        DistMatrix<F,STAR,MC  > z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC> z1_STAR_MC(g), z2_STAR_MC(g);

        z_STAR_MC.AlignWith( U );
        Zeros( z_STAR_MC, 1, m );

        Int k=0;
        while( k < m )
        {
            const Int nbProp = Min(bsize,m-k);
            const bool in2x2 =
                ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
            const Int nb = ( in2x2 ? nbProp+1 : nbProp );

            LockedViewRange( U11, U, k, k,    k+nb, k+nb );
            LockedViewRange( U12, U, k, k+nb, k+nb, m    );

            ViewRange( x1, x, 0, k, 1, k+nb );

            ViewRange( z1_STAR_MC, z_STAR_MC, 0, k,    1, k+nb );
            ViewRange( z2_STAR_MC, z_STAR_MC, 0, k+nb, 1, m    );

            if( k != 0 )
            {
                z1_MR_MC.ColSumScatterFrom( z1_STAR_MC );
                z1.AlignWith( x1 );
                z1 = z1_MR_MC;
                Axpy( F(1), z1, x1 );
            }

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            QuasiTrsvUT
            ( TRANSPOSE, U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(),
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_STAR_MR.AlignWith( U12 );
            x1_STAR_MR = x1_STAR_STAR;
            LocalGemv( TRANSPOSE, F(-1), U12, x1_STAR_MR, F(1), z2_STAR_MC );

            k += nb;
        }
    }
    if( conjugate )
        Conjugate( x );
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_QUASITRSV_UT_HPP
