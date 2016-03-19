/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_FORSYTHE_HPP
#define ELEM_FORSYTHE_HPP

#include "./Jordan.hpp"

namespace elem {

template<typename T> 
inline void
MakeForsythe( Matrix<T>& J, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("MakeForsythe"))
    MakeJordan( J, lambda );
    const Int m = J.Height();
    const Int n = J.Width();
    if( m > 0 && n > 0 )
        J.Set( m-1, 0, alpha );
}

template<typename T,Dist U,Dist V>
inline void
MakeForsythe( DistMatrix<T,U,V>& J, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("MakeForsythe"))
    MakeJordan( J, lambda );
    const Int m = J.Height();
    const Int n = J.Width();
    if( m > 0 && n > 0 )
        J.Set( m-1, 0, alpha );
}

template<typename T,Dist U,Dist V>
inline void
MakeForsythe( BlockDistMatrix<T,U,V>& J, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("MakeForsythe"))
    MakeJordan( J, lambda );
    const Int m = J.Height();
    const Int n = J.Width();
    if( m > 0 && n > 0 )
        J.Set( m-1, 0, alpha );
}

template<typename T>
inline void
Forsythe( Matrix<T>& J, Int n, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Forsythe"))
    J.Resize( n, n );
    MakeForsythe( J, alpha, lambda );
}

template<typename T,Dist U,Dist V>
inline void
Forsythe( DistMatrix<T,U,V>& J, Int n, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Forsythe"))
    J.Resize( n, n );
    MakeForsythe( J, alpha, lambda );
}

template<typename T,Dist U,Dist V>
inline void
Forsythe( BlockDistMatrix<T,U,V>& J, Int n, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Forsythe"))
    J.Resize( n, n );
    MakeForsythe( J, alpha, lambda );
}

} // namespace elem

#endif // ifndef ELEM_FORSYTHE_HPP
