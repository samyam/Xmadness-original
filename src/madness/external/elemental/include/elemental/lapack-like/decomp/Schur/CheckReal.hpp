/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SCHUR_CHECkREAL_HPP
#define ELEM_SCHUR_CHECkREAL_HPP

namespace elem {
namespace schur {

template<typename Real>
inline void
CheckRealSchur( const Matrix<Real>& U, bool standardForm=false )
{
    DEBUG_ONLY(CallStackEntry cse("CheckRealSchur")) 
    const Int n = U.Height();

    auto uMain = U.GetDiagonal();
    auto uSub = U.GetDiagonal( -1 );
    if( standardForm )
    {
        auto uSup = U.GetDiagonal( +1 );
        for( Int j=0; j<n-1; ++j )
        {
            const Real thisDiag = uMain.Get(j,  0);
            const Real nextDiag = uMain.Get(j+1,0);
            const Real thisSub = uSub.Get(j,0);
            const Real thisSup = uSup.Get(j,0);
            if( uSub.Get(j,0) != Real(0) && thisDiag != nextDiag ) 
                LogicError
                ("Diagonal of 2x2 block was not constant: ",thisDiag," and ",
                 nextDiag);
            if( thisSub*thisSup >= 0 )
                LogicError("b*c >= 0: b=",thisSup," and c=",thisSub);
        }
    }

    if( n < 3 )
        return;
    for( Int j=0; j<n-2; ++j )
    {
        const Real thisSub = uSub.Get(j,  0);
        const Real nextSub = uSub.Get(j+1,0);
        if( thisSub != Real(0) && nextSub != Real(0) )
            LogicError
            ("Quasi-triangular assumption broken at j=",j,
             ": subdiagonals were ",thisSub," and ",nextSub);
    }
}

template<typename Real>
inline void
CheckRealSchur( const Matrix<Complex<Real>>& U, bool standardForm=false )
{
    DEBUG_ONLY(CallStackEntry cse("CheckRealSchur")) 
    LogicError("ChceckRealSchur called for complex matrix");
}

template<typename Real>
inline void
CheckRealSchur( const DistMatrix<Real>& U, bool standardForm=false )
{
    DEBUG_ONLY(CallStackEntry cse("CheckRealSchur")) 
    const Int n = U.Height();

    auto uMain = U.GetDiagonal();
    auto uSub = U.GetDiagonal( -1 );
    DistMatrix<Real,STAR,STAR> uMain_STAR_STAR( uMain ),
                               uSub_STAR_STAR( uSub );
    if( standardForm )
    {
        auto uSup = U.GetDiagonal( +1 );
        DistMatrix<Real,STAR,STAR> uSup_STAR_STAR( uSup );
        for( Int j=0; j<n-1; ++j )
        {
            const Real thisDiag = uMain_STAR_STAR.Get(j,  0);
            const Real nextDiag = uMain_STAR_STAR.Get(j+1,0);
            const Real thisSub = uSub_STAR_STAR.Get(j,0);
            const Real thisSup = uSup_STAR_STAR.Get(j,0);
            if( thisSub != Real(0) && thisDiag != nextDiag ) 
                LogicError
                ("Diagonal of 2x2 block was not constant: ",thisDiag," and ",
                 nextDiag);
            if( thisSub*thisSup >= 0 )
                LogicError("b*c >= 0: b=",thisSup," and c=",thisSub);
        }
    }

    if( n < 3 )
        return;
    for( Int j=0; j<n-2; ++j )
    {
        const Real thisSub = uSub_STAR_STAR.Get(j,  0);
        const Real nextSub = uSub_STAR_STAR.Get(j+1,0);
        if( thisSub != Real(0) && nextSub != Real(0) )
            LogicError
            ("Quasi-triangular assumption broken at j=",j,
             ": subdiagonals were ",thisSub," and ",nextSub);
    }
}

template<typename Real>
inline void
CheckRealSchur( const DistMatrix<Complex<Real>>& U, bool standardForm=false )
{
    DEBUG_ONLY(CallStackEntry cse("CheckRealSchur")) 
    LogicError("ChceckRealSchur called for complex matrix");
}

} // namespace schur
} // namespace elem

#endif // ifndef ELEM_SCHUR_CHECkREAL_HPP
