/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_FOURIER_INC
#include ELEM_UNIFORM_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--n","size of matrix",20);
        const Int numRows = Input("--numRows","num rows of submatrix",5);
        const Int numCols = Input("--numCols","num cols of submatrix",5);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        if( numRows > n || numCols > n )
            LogicError("Submatrix too large");

        DistMatrix<Complex<double>> A;
        Fourier( A, n );
        if( display )
            Display( A, "Fourier Matrix" );
        if( print )
            Print( A, "Fourier matrix:" );

        // Get a consistent set of row and column indices (duplication is okay)
        std::vector<Int> rowInds(numRows), colInds(numCols);
        if( mpi::WorldRank() == 0 )
        {
            for( Int j=0; j<numRows; ++j )
                rowInds[j] = SampleUniform<Int>(0,n);
            for( Int j=0; j<numCols; ++j )
                colInds[j] = SampleUniform<Int>(0,n);
        }
        mpi::Broadcast( rowInds.data(), numRows, 0, mpi::COMM_WORLD );
        mpi::Broadcast( colInds.data(), numCols, 0, mpi::COMM_WORLD );
        if( mpi::WorldRank() == 0 && print )
        {
            std::cout << "rowInds: \n";
            for( Int j=0; j<numRows; ++j )
                std::cout << rowInds[j] << "\n";
            std::cout << "\n";
            std::cout << "colInds: \n";
            for( Int j=0; j<numCols; ++j )
                std::cout << colInds[j] << "\n";
            std::cout << std::endl;
        }
        
        auto ASub = A.GetSubmatrix( rowInds, colInds );
        if( display )
            Display( ASub, "ASub" );
        if( print )
            Print( ASub, "ASub" );
      
        MakeUniform( ASub );
        if( display )
            Display( ASub, "Scrambled ASub" );
        if( print )
            Print( ASub, "Scrambled ASub" );
        A.SetSubmatrix( rowInds, colInds, ASub );
       
        if( display )
            Display( A, "Modified Fourier matrix" );
        if( print )
            Print( A, "Modified Fourier matrix" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
