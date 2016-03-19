/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRDTRMM_UNBLOCKED_HPP
#define ELEM_TRDTRMM_UNBLOCKED_HPP

#include ELEM_SYMMETRIC2X2INV_INC
#include ELEM_TRR_INC

namespace elem {
namespace internal {

template<typename F>
inline void
TrdtrmmLUnblocked( Matrix<F>& L, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrdtrmmLUnblocked");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Int n = L.Height();
    const Int ldim = L.LDim();

    Matrix<F> s10;

    for( Int k=0; k<n; ++k )
    {
        auto L00 = ViewRange( L, 0, 0, k,   k );
        auto l10 = ViewRange( L, k, 0, k+1, k );

        // S10 := L10
        s10 = l10;

        // L10 := L10 / delta11
        const F deltaInv = F(1)/L.Get(k,k);
        Scale( deltaInv, l10 );

        // L00 := L00 + l10' s10
        const F* l10Buf = l10.LockedBuffer();
        if( conjugate )
        {
            for( Int j=0; j<k; ++j )
            {
                F* L00Col = L00.Buffer(0,j);
                const F gamma = s10.Get(0,j);
                for( Int i=j; i<k; ++i )
                    L00Col[i] += Conj(l10Buf[i*ldim])*gamma;
            }
        }
        else
        {
            for( Int j=0; j<k; ++j )
            {
                F* L00Col = L00.Buffer(0,j);
                const F gamma = s10.Get(0,j);
                for( Int i=j; i<k; ++i )
                    L00Col[i] += l10Buf[i*ldim]*gamma;
            }
        }

        // L11 := 1 / delta11
        L.Set( k, k, deltaInv );
    }
}

template<typename F>
inline void
TrdtrmmLUnblocked( Matrix<F>& L, const Matrix<F>& dSub, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrdtrmmLUnblocked");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Int n = L.Height();
    const Int ldim = L.LDim();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> s10, S10, D11(2,2);

    Int k=0;
    while( k < n )
    {
        const Int nb = ( k<n-1 && dSub.Get(k,0) != F(0) ? 2 : 1 );

        if( nb == 1 )
        {
            auto L00 = ViewRange( L, 0, 0, k,    k );
            auto l10 = ViewRange( L, k, 0, k+nb, k );

            // S10 := L10
            s10 = l10;

            // L10 := L10 / delta11
            const F deltaInv = F(1)/L.Get(k,k);
            Scale( deltaInv, l10 );

            // L00 := L00 + l10' s10
            // TODO: Extend Trr for this case and then switch
            const F* l10Buf = l10.LockedBuffer();
            if( conjugate )
            {
                for( Int j=0; j<k; ++j )
                {
                    F* L00Col = L00.Buffer(0,j);
                    const F gamma = s10.Get(0,j);
                    for( Int i=j; i<k; ++i )
                        L00Col[i] += Conj(l10Buf[i*ldim])*gamma;
                }
            }
            else
            {
                for( Int j=0; j<k; ++j )
                {
                    F* L00Col = L00.Buffer(0,j);
                    const F gamma = s10.Get(0,j);
                    for( Int i=j; i<k; ++i )
                        L00Col[i] += l10Buf[i*ldim]*gamma;
                }
            }

            // lambda11 := 1 / delta11
            L.Set( k, k, deltaInv );
        }
        else
        {
            auto L00 = ViewRange( L, 0, 0, k,    k    );
            auto L10 = ViewRange( L, k, 0, k+nb, k    );
            auto L11 = ViewRange( L, k, k, k+nb, k+nb );

            // S10 := L10
            S10 = L10;

            // L10 := inv(D11) L10 
            D11.Set( 0, 0, L11.Get(0,0) );
            D11.Set( 1, 1, L11.Get(1,1) );
            D11.Set( 1, 0, dSub.Get(k,0) );
            Symmetric2x2Solve( LEFT, LOWER, D11, L10, conjugate );

            // L00 := L00 + L10' S10
            Trrk( LOWER, orientation, NORMAL, F(1), L10, S10, F(1), L00 );

            // L11 := inv(D11)
            Symmetric2x2Inv( LOWER, D11, conjugate );
            L11.Set( 0, 0, D11.Get(0,0) );
            L11.Set( 1, 0, D11.Get(1,0) );
            L11.Set( 1, 1, D11.Get(1,1) );
        }

        k += nb;
    }
}

template<typename F>
inline void
TrdtrmmUUnblocked( Matrix<F>& U, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrdtrmmUUnblocked");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
    )
    const Int n = U.Height();

    Matrix<F> s01;

    for( Int k=0; k<n; ++k )
    {
        auto U00 = ViewRange( U, 0, 0, k, k   );
        auto u01 = ViewRange( U, 0, k, k, k+1 );

        s01 = u01;

        // u01 := u01 / delta11
        const F deltaInv = F(1)/U.Get(k,k);
        Scale( deltaInv, u01 );

        // U00 := U00 + s01 u01'
        Trr( UPPER, F(1), s01, u01, U00, conjugate );

        // lambda11 := 1 / delta11
        U.Set( k, k, deltaInv );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRDTRMM_UNBLOCKED_HPP
