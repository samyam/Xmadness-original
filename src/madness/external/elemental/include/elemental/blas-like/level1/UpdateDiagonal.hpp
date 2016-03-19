/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_UPDATEDIAGONAL_HPP
#define ELEM_UPDATEDIAGONAL_HPP

// This is essentially equivalent to SetDiagonal, but with s/Set/Update/g.

namespace elem {

template<typename T,typename S>
inline void
UpdateDiagonal( Matrix<T>& A, S alpha )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    ELEM_PARALLEL_FOR
    for( Int j=0; j<Min(height,width); ++j )
        A.Update(j,j,alpha);
}

template<typename T,typename S>
inline void
UpdateDiagonal( Matrix<T>& A, S alpha, Int offset, LeftOrRight side=LEFT )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( side == LEFT )
    {
        for( Int j=0; j<width; ++j )
        {
            const Int i = j-offset;
            if( i >= 0 && i < height )
                A.Update(i,j,alpha);
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            const Int i = j-offset+height-width;
            if( i >= 0 && i < height )
                A.Update(i,j,alpha);
        }
    }
}

template<typename T,typename S,Dist U,Dist V>
inline void
UpdateDiagonal( DistMatrix<T,U,V>& A, S alpha )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateDiagonal"))
    const Int height = A.Height();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        if( j < height && A.IsLocalRow(j) )
        {
            const Int iLoc = A.LocalRow(j);
            A.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,typename S,Dist U,Dist V>
inline void
UpdateDiagonal
( DistMatrix<T,U,V>& A, S alpha, Int offset, LeftOrRight side=LEFT )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int localWidth = A.LocalWidth();
    if( side == LEFT )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int i = j-offset;
            if( i >= 0 && i < height && A.IsLocalRow(i) )
            {
                const Int iLoc = A.LocalRow(i);
                A.UpdateLocal( iLoc, jLoc, alpha );
            }
        }
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int i = j-offset+height-width;
            if( i >= 0 && i < height && A.IsLocalRow(i) )
            {
                const Int iLoc = A.LocalRow(i);
                A.UpdateLocal( iLoc, jLoc, alpha );
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_UPDATEDIAGONAL_HPP
