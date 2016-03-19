/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IMPORTS_FLAME_HPP
#define ELEM_IMPORTS_FLAME_HPP

#ifdef ELEM_HAVE_FLA_BSVD
namespace elem {

void FlaBidiagSVD
( int k, int mU, int mV, double* d, double* e, 
  double* U, int ldu, double* V, int ldv, 
  int numAccum=32, int maxNumIts=30, int bAlg=512 );

void FlaBidiagSVD
( int k, int mU, int mV, double* d, double* e, 
  Complex<double>* U, int ldu, Complex<double>* V, int ldv, 
  int numAccum=32, int maxNumIts=30, int bAlg=512 );

} // namespace elem
#endif // ifdef ELEM_HAVE_FLA_BSVD

#endif // ifndef ELEM_IMPORTS_FLAME_HPP
