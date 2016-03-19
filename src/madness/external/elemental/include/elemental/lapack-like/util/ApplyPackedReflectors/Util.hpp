/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_APPLYPACKEDREFLECTORS_UTIL_HPP
#define ELEM_APPLYPACKEDREFLECTORS_UTIL_HPP

namespace elem {

template<typename F> 
inline void
FixDiagonal
( Conjugation conjugation, const Matrix<F>& t, Matrix<F>& SInv )
{
    DEBUG_ONLY(CallStackEntry cse("FixDiagonal"))
    for( Int j=0; j<SInv.Height(); ++j )
    {
        const F value = t.Get(j,0);
        if( conjugation == CONJUGATED )
            SInv.Set(j,j,F(1)/Conj(value));
        else
            SInv.Set(j,j,F(1)/value);
    }
}

template<typename F> 
inline void
FixDiagonal
( Conjugation conjugation,
  const DistMatrix<F,STAR,STAR>& t,
        DistMatrix<F,STAR,STAR>& SInv )
{
    DEBUG_ONLY(CallStackEntry cse("FixDiagonal"))
    for( Int j=0; j<SInv.Height(); ++j )
    {
        const F value = t.GetLocal(j,0);
        if( conjugation == CONJUGATED )
            SInv.SetLocal( j, j, F(1)/Conj(value) );
        else
            SInv.SetLocal( j, j, F(1)/value );
    }
}

} // namespace elem

#endif // ifndef ELEM_APPLYPACKEDREFLECTORS_UTIL_HPP
