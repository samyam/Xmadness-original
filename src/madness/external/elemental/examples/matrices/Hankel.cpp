/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HANKEL_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--height","height of matrix",10);
        const Int n = Input("--width","width of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        const Int length = m+n-1;
        std::vector<double> a( length );
        for( Int j=0; j<length; ++j )
            a[j] = j;

        DistMatrix<double> H;
        Hankel( H, m, n, a );
        if( display )
            Display( H, "Hankel" );
        if( print )
            Print( H, "Hankel matrix:" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
