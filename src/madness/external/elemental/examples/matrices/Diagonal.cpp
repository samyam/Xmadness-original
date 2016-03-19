/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_DIAGONAL_INC
#include ELEM_ZEROS_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        std::vector<double> d( n );
        for( Int j=0; j<n; ++j )
            d[j] = j;

        DistMatrix<double> D;
        Diagonal( D, d );
        if( display )
        {
            Display( D, "Diagonal matrix" );
#ifdef ELEM_HAVE_QT5
            Spy( D, "Diagonal spy plot" );
#endif
        }
        if( print )
            Print( D, "D:" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
