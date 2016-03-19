/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--height","height of matrix",10);
        const Int n = Input("--width","width of matrix",10);
        const std::string filename = 
            Input("--filename","filename",std::string(""));
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        if( filename == "" )
            LogicError("Please specify a filename to read");

        DistMatrix<double> A(m,n);
        Read( A, filename );
        if( display )
            Display( A, "A (distributed read)" );
        if( print )
            Print( A, "A (distributed read)" );
        Read( A, filename, AUTO, true );
        if( display )
            Display( A, "A (sequential read)" );
        if( print )
            Print( A, "A (sequential read)" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
