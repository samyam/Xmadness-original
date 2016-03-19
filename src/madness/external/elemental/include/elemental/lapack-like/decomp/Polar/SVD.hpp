/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_POLAR_SVD_HPP
#define ELEM_POLAR_SVD_HPP

#include ELEM_GEMM_INC

#include ELEM_SVD_INC

#include ELEM_HERMITIANFROMEVD_INC
#include ELEM_ZEROS_INC

namespace elem {
namespace polar {

// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.

template<typename F>
inline void
SVD( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("polar::SVD"))
    // Get the SVD of A
    typedef Base<F> Real;
    Matrix<Real> s;
    Matrix<F> U, V;
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );
}

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<F>& P )
{
    DEBUG_ONLY(CallStackEntry cse("polar::SVD"))
    // Get the SVD of A
    typedef Base<F> Real;
    Matrix<Real> s;
    Matrix<F> U, V;
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

template<typename F>
inline void
SVD( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("polar::SVD"))
    // Get the SVD of A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );
}

template<typename F>
inline void
SVD( DistMatrix<F>& A, DistMatrix<F>& P )
{
    DEBUG_ONLY(CallStackEntry cse("polar::SVD"))
    // Get the SVD of A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

} // namespace polar
} // namespace elem

#endif // ifndef ELEM_POLAR_SVD_HPP
