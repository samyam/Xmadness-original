/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRTRMM_UNBLOCKED_HPP
#define ELEM_TRTRMM_UNBLOCKED_HPP

namespace elem {
namespace internal {

template<typename T>
inline void
TrtrmmLUnblocked( Matrix<T>& L, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrtrmmLUnblocked");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Int n = L.Height();

    T* LBuffer = L.Buffer();
    const Int ldim = L.LDim();
    for( Int j=0; j<n; ++j )
    {
        T* ELEM_RESTRICT l10 = &LBuffer[j];
        if( conjugate )
        {
            // L00 := L00 + l10^H l10
            for( Int k=0; k<j; ++k )
            {
                const T gamma = l10[k*ldim];
                T* ELEM_RESTRICT L00Col = &LBuffer[k*ldim];
                for( Int i=k; i<j; ++i )
                    L00Col[i] += Conj(l10[i*ldim])*gamma;
            }
        }
        else
        {
            // L00 := L00 + l10^T l10
            for( Int k=0; k<j; ++k )
            {
                const T gamma = l10[k*ldim];
                T* ELEM_RESTRICT L00Col = &LBuffer[k*ldim];
                for( Int i=k; i<j; ++i )
                    L00Col[i] += l10[i*ldim]*gamma;
            }
        }

        // l10 := l10 lambda11
        const T lambda11 = LBuffer[j+j*ldim];
        for( Int k=0; k<j; ++k )
            l10[k*ldim] *= lambda11;

        // lambda11 := lambda11^2 or |lambda11|^2
        if( conjugate )
            LBuffer[j+j*ldim] = lambda11*Conj(lambda11);
        else
            LBuffer[j+j*ldim] = lambda11*lambda11;
    }
}

template<typename T>
inline void
TrtrmmUUnblocked( Matrix<T>& U, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrtrmmUUnblocked");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
    )
    const Int n = U.Height();

    T* UBuffer = U.Buffer();
    const Int ldim = U.LDim();
    for( Int j=0; j<n; ++j )
    {
        T* ELEM_RESTRICT u01 = &UBuffer[j*ldim];
        if( conjugate )
        {
            // U00 := U00 + u01 u01^H
            for( Int k=0; k<j; ++k )
            {
                const T gamma = Conj(u01[k]);
                T* ELEM_RESTRICT U00Col = &UBuffer[k*ldim];
                for( Int i=0; i<=k; ++i )
                    U00Col[i] += u01[i]*gamma;
            }
        }
        else
        {
            // U00 := U00 + u01 u01^T
            for( Int k=0; k<j; ++k )
            {
                const T gamma = u01[k];
                T* ELEM_RESTRICT U00Col = &UBuffer[k*ldim];
                for( Int i=0; i<=k; ++i )
                    U00Col[i] += u01[i]*gamma;
            }
        }

        // u01 := u01 upsilon11
        const T upsilon11 = UBuffer[j+j*ldim];
        for( Int k=0; k<j; ++k )
            u01[k] *= upsilon11;

        // upsilon11 := upsilon11^2 or |upsilon11|^2
        if( conjugate )
            UBuffer[j+j*ldim] = upsilon11*Conj(upsilon11);
        else
            UBuffer[j+j*ldim] = upsilon11*upsilon11;
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRTRMM_UNBLOCKED_HPP
