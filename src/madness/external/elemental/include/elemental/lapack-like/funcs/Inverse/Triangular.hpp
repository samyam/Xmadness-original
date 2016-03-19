/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_INVERSE_TRIANGULAR_HPP
#define ELEM_INVERSE_TRIANGULAR_HPP

namespace elem {
template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );
} // namespace elem

#include "./Triangular/LVar3.hpp"
#include "./Triangular/UVar3.hpp"

namespace elem {
namespace triang_inv {

template<typename F>
inline void
Var3( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A  )
{
    DEBUG_ONLY(CallStackEntry cse("triang_inv::Var3"))
    if( uplo == LOWER )
        LVar3( diag, A );
    else
        UVar3( diag, A );
}

template<typename F>
inline void
Var3( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  )
{
    DEBUG_ONLY(CallStackEntry cse("triang_inv::Var3"))
    if( uplo == LOWER )
        LVar3( diag, A );
    else
        UVar3( diag, A );
}

} // namespace triang_inv

template<typename F>
inline void
TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularInverse"))
    triang_inv::Var3( uplo, diag, A );
}

template<typename F>
inline void
TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularInverse"))
    triang_inv::Var3( uplo, diag, A );
}

template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTriangularInverse"))
    TriangularInverse( uplo, diag, A.Matrix() );
}

} // namespace elem

#endif // ifndef ELEM_INVERSE_TRIANGULAR_HPP
