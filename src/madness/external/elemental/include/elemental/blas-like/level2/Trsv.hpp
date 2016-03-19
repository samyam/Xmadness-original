/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSV_HPP
#define ELEM_TRSV_HPP

#include "./Trsv/LN.hpp"
#include "./Trsv/LT.hpp"
#include "./Trsv/UN.hpp"
#include "./Trsv/UT.hpp"

namespace elem {

template<typename F>
inline void
Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<F>& A, Matrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trsv");
        if( x.Height() != 1 && x.Width() != 1 )
            LogicError("x must be a vector");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        if( xLength != A.Height() )
            LogicError("x must conform with A");
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    const Int m = A.Height();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
}

template<typename F>
inline void
Trsv
( UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  const DistMatrix<F>& A,
        DistMatrix<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Trsv"))
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrsvLN( diag, A, x );
        else
            internal::TrsvLT( orientation, diag, A, x );
    }
    else
    {
        if( orientation == NORMAL )
            internal::TrsvUN( diag, A, x );
        else
            internal::TrsvUT( orientation, diag, A, x );
    }
}

} // namespace elem

#endif // ifndef ELEM_TRSV_HPP
