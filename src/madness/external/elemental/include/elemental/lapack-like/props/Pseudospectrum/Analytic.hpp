/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_ANALYTIC_HPP
#define ELEM_PSEUDOSPECTRUM_ANALYTIC_HPP

#include "./Util.hpp"

namespace elem {
namespace pspec {

template<typename Real>
inline void
Analytic
( const Matrix<Complex<Real>>& w, 
  const Matrix<Complex<Real>>& shifts, 
        Matrix<Real         >& invNorms,
        SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Analytic"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = w.Height();
    const Int numShifts = shifts.Height();
    const Real normCap = NormCap<Real>();

    Zeros( invNorms, numShifts, 1 );
    if( n == 0 )
        return;

    for( Int j=0; j<numShifts; ++j )
    {
        const C shift = shifts.Get(j,0);
        Real minDist = Abs(shift-w.Get(0,0));
        for( Int k=1; k<n; ++k )
        {
            const Real dist = Abs(shift-w.Get(k,0));
            minDist = Min(dist,minDist);
        }
        Real alpha = Real(1)/minDist;
        if( std::isnan(alpha) || alpha >= normCap )
            alpha = normCap;
        invNorms.Set( j, 0, alpha );
    }
    
    snapCtrl.itCounts = false;
    Matrix<Int> itCounts;
    FinalSnapshot( invNorms, itCounts, snapCtrl );
}

template<typename Real,Dist colDist,Dist rowDist>
inline void
Analytic
( const DistMatrix<Complex<Real>,colDist,rowDist>& w, 
  const DistMatrix<Complex<Real>,VR,     STAR   >& shifts, 
        DistMatrix<Real,         VR,     STAR   >& invNorms,
        SnapshotCtrl& snapCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Analytic"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = w.Height();
    const Int numShifts = shifts.Height();
    const Real normCap = NormCap<Real>();

    Zeros( invNorms, numShifts, 1 );
    if( n == 0 )
        return;

    DistMatrix<C,STAR,STAR> w_STAR_STAR( w );

    const Int numLocShifts = shifts.LocalHeight();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const C shift = shifts.GetLocal(jLoc,0);
        Real minDist = Abs(shift-w_STAR_STAR.GetLocal(0,0));
        for( Int k=1; k<n; ++k )
        {
            const Real dist = Abs(shift-w_STAR_STAR.GetLocal(k,0));
            minDist = Min(dist,minDist);
        }
        Real alpha = Real(1)/minDist;
        if( std::isnan(alpha) || alpha >= normCap )
            alpha = normCap;
        invNorms.SetLocal( jLoc, 0, alpha );
    }

    snapCtrl.itCounts = false;
    DistMatrix<Int,VR,STAR> itCounts(w.Grid());
    FinalSnapshot( invNorms, itCounts, snapCtrl );
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_ANALYTIC_HPP
