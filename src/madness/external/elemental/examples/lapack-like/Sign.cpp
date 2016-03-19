/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_SIGN_INC
#include ELEM_UNIFORM_INC
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const SignScaling scaling = 
            static_cast<SignScaling>(Input("--scaling","scaling strategy",0));
        const Int maxIts = Input("--maxIts","max number of iter's",100);
        const double tol = Input("--tol","convergence tolerance",1e-6);
        const bool progress = Input("--progress","print sign progress?",true);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );

        SignCtrl<Real> signCtrl;
        signCtrl.maxIts = maxIts;
        signCtrl.tol = tol;
        signCtrl.progress = progress;
        signCtrl.scaling = scaling;

        // Compute sgn(A)
        Sign( A, signCtrl );
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
