/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HELMHOLTZPML_HPP
#define ELEM_HELMHOLTZPML_HPP

#include "./Zeros.hpp"

namespace elem {

namespace pml {

template<typename Real>
Complex<Real> 
Profile( Real x, Real w, Real pmlExp, Real sigma, Real k )
{
    DEBUG_ONLY(
        CallStackEntry cse("pml::Profile");
        if( x < 0 || x > w )
            LogicError("Evaluation point not in PML interval");
    )
    const Real realPart(1);
    const Real arg = x/w;
    const Real imagPart = (sigma/w)*Pow(arg,pmlExp)/k;
    return Complex<Real>(realPart,imagPart);
}

template<typename Real>
Complex<Real>
sInv( Int j, Int n, Int numPmlPoints, Real h, Real pmlExp, Real sigma, Real k )
{
    if( j < numPmlPoints-1 )
        return Profile
               ( (numPmlPoints-1-j)*h, numPmlPoints*h, pmlExp, sigma, k );
    else if( j > n-numPmlPoints )
        return Profile
               ( (j-(n-numPmlPoints))*h, numPmlPoints*h, pmlExp, sigma, k );
    else
        return Complex<Real>(1,0);
}

} // namespace pml

// 1D Helmholtz with PML
template<typename Real> 
inline void
HelmholtzPML
( Matrix<Complex<Real>>& H, Int n, Complex<Real> omega, 
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    using namespace pml;
    typedef Complex<Real> C;
    Zeros( H, n, n );

    const Real k = RealPart(omega) / (2*M_PI);
    const Real h = Real(1)/(n+1);
    const Real hSquared = h*h;
    for( Int i=0; i<n; ++i )
    {
        const Int x = i;

        const C sInvL = sInv( x-1, n, numPmlPoints, h, pmlExp, sigma, k );
        const C sInvM = sInv( x,   n, numPmlPoints, h, pmlExp, sigma, k );
        const C sInvR = sInv( x+1, n, numPmlPoints, h, pmlExp, sigma, k );

        // This is a bit silly in 1D, but it keeps the structure from 2D/3D
        const C xTop = Real(1);
        const C xTempL = xTop/sInvL;
        const C xTempM = xTop/sInvM;
        const C xTempR = xTop/sInvR;
        const C xTermL = (xTempL+xTempM) / (2*hSquared);
        const C xTermR = (xTempM+xTempR) / (2*hSquared);

        const C mainTerm = (xTermL+xTermR) - omega*omega*sInvM;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -xTermL );
        if( x != n-1 )
            H.Set( i, i+1, -xTermR );
    }
}

// 1D Helmholtz with PML
template<typename Real> 
inline void
HelmholtzPML
( DistMatrix<Complex<Real>>& H, Int n, Complex<Real> omega, 
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    using namespace pml;
    typedef Complex<Real> C;
    Zeros( H, n, n );

    const Real k = RealPart(omega) / (2*M_PI);
    const Real h = Real(1)/(n+1);
    const Real hSquared = h*h;

    const Int localHeight = H.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i;

        const C sInvL = sInv( x-1, n, numPmlPoints, h, pmlExp, sigma, k );
        const C sInvM = sInv( x,   n, numPmlPoints, h, pmlExp, sigma, k );
        const C sInvR = sInv( x+1, n, numPmlPoints, h, pmlExp, sigma, k );

        // This is a bit silly in 1D, but it keeps the structure from 2D/3D
        const C xTop = Real(1);
        const C xTempL = xTop/sInvL;
        const C xTempM = xTop/sInvM;
        const C xTempR = xTop/sInvR;
        const C xTermL = (xTempL+xTempM) / (2*hSquared);
        const C xTermR = (xTempM+xTempR) / (2*hSquared);

        const C mainTerm = (xTermL+xTermR) - omega*omega*sInvM;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -xTermL );
        if( x != n-1 )
            H.Set( i, i+1, -xTermR );
    }
}

// 2D Helmholtz with PML
template<typename Real> 
inline void
HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Int ny, Complex<Real> omega, 
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    using namespace pml;
    typedef Complex<Real> C;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const Real k = RealPart(omega) / (2*M_PI);
    const Real hx = Real(1)/(nx+1);
    const Real hy = Real(1)/(ny+1);
    const Real hxSquared = hx*hx;
    const Real hySquared = hy*hy;
    for( Int i=0; i<n; ++i )
    {
        const Int x = i % nx;
        const Int y = i / nx; 

        const C sxInvL = sInv( x-1, nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvM = sInv( x,   nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvR = sInv( x+1, nx, numPmlPoints, hx, pmlExp, sigma, k );

        const C syInvL = sInv( y-1, ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvM = sInv( y,   ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvR = sInv( y+1, ny, numPmlPoints, hy, pmlExp, sigma, k );

        const C xTop = syInvM;
        const C xTempL = xTop/sxInvL;
        const C xTempM = xTop/sxInvM;
        const C xTempR = xTop/sxInvR;
        const C xTermL = (xTempL+xTempM) / (2*hxSquared);
        const C xTermR = (xTempM+xTempR) / (2*hxSquared);

        const C yTop = sxInvM;
        const C yTempL = yTop/syInvL;
        const C yTempM = yTop/syInvM;
        const C yTempR = yTop/syInvR;
        const C yTermL = (yTempL+yTempM) / (2*hySquared);
        const C yTermR = (yTempM+yTempR) / (2*hySquared);

        const C mainTerm = (xTermL+xTermR+yTermL+yTermR) - 
                           omega*omega*sxInvM*syInvM;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -xTermL );
        if( x != nx-1 )
            H.Set( i, i+1, -xTermR );
        if( y != 0 )
            H.Set( i, i-nx, -yTermL );
        if( y != ny-1 )
            H.Set( i, i+nx, -yTermR );
    }
}

// 2D Helmholtz with PML
template<typename Real> 
inline void
HelmholtzPML
( DistMatrix<Complex<Real>>& H, Int nx, Int ny, Complex<Real> omega, 
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    using namespace pml;
    typedef Complex<Real> C;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const Real k = RealPart(omega) / (2*M_PI);
    const Real hx = Real(1)/(nx+1);
    const Real hy = Real(1)/(ny+1);
    const Real hxSquared = hx*hx;
    const Real hySquared = hy*hy;

    const Int localHeight = H.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i % nx;
        const Int y = i / nx; 

        const C sxInvL = sInv( x-1, nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvM = sInv( x,   nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvR = sInv( x+1, nx, numPmlPoints, hx, pmlExp, sigma, k );

        const C syInvL = sInv( y-1, ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvM = sInv( y,   ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvR = sInv( y+1, ny, numPmlPoints, hy, pmlExp, sigma, k );

        const C xTop = syInvM;
        const C xTempL = xTop/sxInvL;
        const C xTempM = xTop/sxInvM;
        const C xTempR = xTop/sxInvR;
        const C xTermL = (xTempL+xTempM) / (2*hxSquared);
        const C xTermR = (xTempM+xTempR) / (2*hxSquared);

        const C yTop = sxInvM;
        const C yTempL = yTop/syInvL;
        const C yTempM = yTop/syInvM;
        const C yTempR = yTop/syInvR;
        const C yTermL = (yTempL+yTempM) / (2*hySquared);
        const C yTermR = (yTempM+yTempR) / (2*hySquared);

        const C mainTerm = (xTermL+xTermR+yTermL+yTermR) - 
                           omega*omega*sxInvM*syInvM;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -xTermL );
        if( x != nx-1 )
            H.Set( i, i+1, -xTermR );
        if( y != 0 )
            H.Set( i, i-nx, -yTermL );
        if( y != ny-1 )
            H.Set( i, i+nx, -yTermR );
    }
}

// 3D Helmholtz with PML
template<typename Real> 
inline void
HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Int ny, Int nz, Complex<Real> omega, 
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    using namespace pml;
    typedef Complex<Real> C;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const Real k = RealPart(omega) / (2*M_PI);
    const Real hx = Real(1)/(nx+1);
    const Real hy = Real(1)/(ny+1);
    const Real hz = Real(1)/(nz+1);
    const Real hxSquared = hx*hx;
    const Real hySquared = hy*hy;
    const Real hzSquared = hz*hz;
    for( Int i=0; i<n; ++i )
    {
        const Int x = i % nx;
        const Int y = (i/nx) % ny; 
        const Int z = i/(nx*ny);

        const C sxInvL = sInv( x-1, nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvM = sInv( x,   nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvR = sInv( x+1, nx, numPmlPoints, hx, pmlExp, sigma, k );

        const C syInvL = sInv( y-1, ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvM = sInv( y,   ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvR = sInv( y+1, ny, numPmlPoints, hy, pmlExp, sigma, k );

        const C szInvL = sInv( z-1, nz, numPmlPoints, hz, pmlExp, sigma, k );
        const C szInvM = sInv( z,   nz, numPmlPoints, hz, pmlExp, sigma, k );
        const C szInvR = sInv( z+1, nz, numPmlPoints, hz, pmlExp, sigma, k );

        const C xTop = syInvM*szInvM;
        const C xTempL = xTop/sxInvL;
        const C xTempM = xTop/sxInvM;
        const C xTempR = xTop/sxInvR;
        const C xTermL = (xTempL+xTempM) / (2*hxSquared);
        const C xTermR = (xTempM+xTempR) / (2*hxSquared);

        const C yTop = sxInvM*szInvM;
        const C yTempL = yTop/syInvL;
        const C yTempM = yTop/syInvM;
        const C yTempR = yTop/syInvR;
        const C yTermL = (yTempL+yTempM) / (2*hySquared);
        const C yTermR = (yTempM+yTempR) / (2*hySquared);

        const C zTop = sxInvM*syInvM;
        const C zTempL = zTop/szInvL;
        const C zTempM = zTop/szInvM;
        const C zTempR = zTop/szInvR;
        const C zTermL = (zTempL+zTempM) / (2*hzSquared);
        const C zTermR = (zTempM+zTempR) / (2*hzSquared);

        const C mainTerm = (xTermL+xTermR+yTermL+yTermR+zTermL+zTermR) - 
                           omega*omega*sxInvM*syInvM*szInvM;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -xTermL );
        if( x != nx-1 )
            H.Set( i, i+1, -xTermR );
        if( y != 0 )
            H.Set( i, i-nx, -yTermL );
        if( y != ny-1 )
            H.Set( i, i+nx, -yTermR );
        if( z != 0 )
            H.Set( i, i-nx*ny, -zTermL );
        if( z != nz-1 )
            H.Set( i, i+nx*ny, -zTermR );
    }
}

// 3D Helmholtz with PML
template<typename Real> 
inline void
HelmholtzPML
( DistMatrix<Complex<Real>>& H, Int nx, Int ny, Int nz, Complex<Real> omega, 
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    using namespace pml;
    typedef Complex<Real> C;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const Real k = RealPart(omega) / (2*M_PI);
    const Real hx = Real(1)/(nx+1);
    const Real hy = Real(1)/(ny+1);
    const Real hz = Real(1)/(nz+1);
    const Real hxSquared = hx*hx;
    const Real hySquared = hy*hy;
    const Real hzSquared = hz*hz;

    const Int localHeight = H.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i % nx;
        const Int y = (i/nx) % ny; 
        const Int z = i/(nx*ny);

        const C sxInvL = sInv( x-1, nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvM = sInv( x,   nx, numPmlPoints, hx, pmlExp, sigma, k );
        const C sxInvR = sInv( x+1, nx, numPmlPoints, hx, pmlExp, sigma, k );

        const C syInvL = sInv( y-1, ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvM = sInv( y,   ny, numPmlPoints, hy, pmlExp, sigma, k );
        const C syInvR = sInv( y+1, ny, numPmlPoints, hy, pmlExp, sigma, k );

        const C szInvL = sInv( z-1, nz, numPmlPoints, hz, pmlExp, sigma, k );
        const C szInvM = sInv( z,   nz, numPmlPoints, hz, pmlExp, sigma, k );
        const C szInvR = sInv( z+1, nz, numPmlPoints, hz, pmlExp, sigma, k );

        const C xTop = syInvM*szInvM;
        const C xTempL = xTop/sxInvL;
        const C xTempM = xTop/sxInvM;
        const C xTempR = xTop/sxInvR;
        const C xTermL = (xTempL+xTempM) / (2*hxSquared);
        const C xTermR = (xTempM+xTempR) / (2*hxSquared);

        const C yTop = sxInvM*szInvM;
        const C yTempL = yTop/syInvL;
        const C yTempM = yTop/syInvM;
        const C yTempR = yTop/syInvR;
        const C yTermL = (yTempL+yTempM) / (2*hySquared);
        const C yTermR = (yTempM+yTempR) / (2*hySquared);

        const C zTop = sxInvM*syInvM;
        const C zTempL = zTop/szInvL;
        const C zTempM = zTop/szInvM;
        const C zTempR = zTop/szInvR;
        const C zTermL = (zTempL+zTempM) / (2*hzSquared);
        const C zTermR = (zTempM+zTempR) / (2*hzSquared);

        const C mainTerm = (xTermL+xTermR+yTermL+yTermR+zTermL+zTermR) - 
                           omega*omega*sxInvM*syInvM*szInvM;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -xTermL );
        if( x != nx-1 )
            H.Set( i, i+1, -xTermR );
        if( y != 0 )
            H.Set( i, i-nx, -yTermL );
        if( y != ny-1 )
            H.Set( i, i+nx, -yTermR );
        if( z != 0 )
            H.Set( i, i-nx*ny, -zTermL );
        if( z != nz-1 )
            H.Set( i, i+nx*ny, -zTermR );
    }
}

} // namespace elem

#endif // ifndef ELEM_HELMHOLTZPML_HPP
