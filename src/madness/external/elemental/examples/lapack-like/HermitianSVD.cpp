/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_SVD_INC
#include ELEM_TRACE_INC
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and 
    // does so if necessary. The full routine is elem::Initialize.
    Initialize( argc, argv );

    // Extract our MPI rank
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    // Surround the Elemental calls with try/catch statements in order to 
    // safely handle any exceptions that were thrown during execution.
    try 
    {
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Create a 2d process grid from a communicator. In our case, it is
        // MPI_COMM_WORLD. There is another constructor that allows you to 
        // specify the grid dimensions, Grid g( comm, r, c ), which creates an 
        // r x c grid.
        Grid g( comm );
    
        // Create an n x n complex distributed matrix.
        // We distribute the matrix using grid 'g'.
        //
        // There are quite a few available constructors, including ones that 
        // allow you to pass in your own local buffer and to specify the 
        // distribution alignments (i.e., which process row and column owns the
        // top-left element)
        const Int n = 6; // choose a small problem size since we will print
        DistMatrix<C> H( n, n, g );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
        // the global matrix is Hermitian. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        //
        const Int localHeight = H.LocalHeight();
        const Int localWidth = H.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            // Our process owns the rows colShift:colStride:n,
            //           and the columns rowShift:rowStride:n
            const Int j = H.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = H.GlobalRow(iLoc);
                H.SetLocal( iLoc, jLoc, C(i+j,i-j) );
            }
        }
        // Alternatively, we could have sequentially filled the matrix with 
        // for( Int j=0; j<A.Width(); ++j )
        //   for( Int i=0; i<A.Height(); ++i )
        //     A.Set( i, j, C(i+j,i-j) );
        //
        // More convenient interfaces are being investigated.
        //
        if( print )
            Print( H, "H" );

        // Print its trace
        const C trace = Trace( H );
        if( commRank == 0 )
            std::cout << "Tr(H) = " << trace << std::endl;

        // Build the singular value decomposition through the Hermitian EVD.
        //
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        //
        DistMatrix<R,VR,STAR> s( g );
        DistMatrix<C> U( g ), V( g );
        HermitianSVD( LOWER, H, s, U, V ); // only use lower half of H
        if( print )
        {
            Print( s, "Singular values of H" );
            Print( U, "Left singular vectors of H" );
            Print( V, "Right singular vectors of H" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}

