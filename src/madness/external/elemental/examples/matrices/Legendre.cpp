/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HERMITIANTRIDIAGEIG_INC
#include ELEM_LEGENDRE_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> J;
        Legendre( J, n );
        if( display )
        {
            Display( J, "Jacobi matrix for Legendre polynomials" );
#ifdef ELEM_HAVE_QT5
            Spy( J, "Spy plot for Jacobi matrix" );
#endif
        }
        if( print )
            Print( J, "Jacobi matrix for Legendre polynomials" );

        // We will compute Gaussian quadrature points and weights over [-1,+1]
        // using the eigenvalue decomposition of the Jacobi matrix for the 
        // Legendre polynomials.
        DistMatrix<double,VR,  STAR> points;
        DistMatrix<double,STAR,VR  > X;
        HermitianTridiagEig
        ( J.GetDiagonal(), J.GetDiagonal(-1), points, X, ASCENDING );
        if( display )
            Display( points, "Quadrature points" );
        if( print )
            Print( points, "points" );
        auto firstRow = View( X, 0, 0, 1, n );
        DistMatrix<double,STAR,STAR> weights = firstRow;
        for( Int j=0; j<n; ++j )
        {
            const double gamma = weights.Get( 0, j );
            weights.Set( 0, j, 2*gamma*gamma );
        }
        if( display )
            Display( weights, "Quadrature weights" );
        if( print )
            Print( weights, "weights" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
