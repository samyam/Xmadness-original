/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LOGBARRIER_HPP
#define ELEM_LOGBARRIER_HPP

#include ELEM_DETERMINANT_INC

namespace elem {

template<typename F>
inline Base<F> 
LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F>
inline Base<F>
LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
inline Base<F>
LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
inline Base<F>
LogBarrier( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

} // namespace elem

#endif // ifndef ELEM_LOGBARRIER_HPP
