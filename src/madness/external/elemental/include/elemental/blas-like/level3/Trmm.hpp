/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRMM_HPP
#define ELEM_TRMM_HPP

namespace elem {

template<typename T,Dist BColDist,Dist BRowDist>
inline void
LocalTrmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T,STAR,STAR>& A,
                 DistMatrix<T,BColDist,BRowDist>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrmm");
        if( (side == LEFT && BColDist != STAR) ||
            (side == RIGHT && BRowDist != STAR) )
            LogicError
            ("Dist of RHS must conform with that of triangle");
    )
    Trmm
    ( side, uplo, orientation, diag, alpha, A.LockedMatrix(), B.Matrix() );
}

} // namespace elem

#include "./Trmm/LLN.hpp"
#include "./Trmm/LLT.hpp"
#include "./Trmm/LUN.hpp"
#include "./Trmm/LUT.hpp"
#include "./Trmm/RLN.hpp"
#include "./Trmm/RLT.hpp"
#include "./Trmm/RUN.hpp"
#include "./Trmm/RUT.hpp"

namespace elem {

template<typename T>
inline void
Trmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trmm");
        if( A.Height() != A.Width() )
            LogicError("Triangular matrix must be square");
        if( side == LEFT )
        {
            if( A.Height() != B.Height() )
                LogicError("Nonconformal Trmm");
        }
        else
        {
            if( A.Height() != B.Width() )
                LogicError("Nonconformal Trmm");
        }
    )
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    blas::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
}

template<typename T>
inline void
Trmm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& A, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("Trmm"))
    Scale( alpha, X );
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            trmm::LLN( diag, A, X );
        else
            trmm::LLT( orientation, diag, A, X );
    }
    else if( side == LEFT )
    {
        if( orientation == NORMAL )
            trmm::LUN( diag, A, X );
        else
            trmm::LUT( orientation, diag, A, X );
    }
    else if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            trmm::RLN( diag, A, X );
        else
            trmm::RLT( orientation, diag, A, X );
    }
    else
    {
        if( orientation == NORMAL )
            trmm::RUN( diag, A, X );
        else
            trmm::RUT( orientation, diag, A, X );
    }
}

} // namespace elem

#endif // ifndef ELEM_TRMM_HPP
