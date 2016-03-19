/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SYR2K_INC
#include ELEM_UNIFORM_INC
using namespace std;
using namespace elem;

template<typename T>
void TestSyr2k
( bool print, UpperOrLower uplo, Orientation orientation,
  Int m, Int k, T alpha, T beta, const Grid& g )
{
    DistMatrix<T> A(g), B(g), C(g);

    if( orientation == NORMAL )
    {
        Uniform( A, m, k );
        Uniform( B, m, k );
    }
    else
    {
        Uniform( A, k, m );
        Uniform( B, k, m );
    }
    Uniform( C, m, m );
    MakeTriangular( uplo, C );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting Syr2k...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Syr2k( uplo, orientation, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2.*double(m)*double(m)*double(k)/(1.e9*runTime);
    const double gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        ostringstream msg;
        if( orientation == NORMAL )
            msg << "C := " << alpha << " A B' + B A'" << beta << " C";
        else
            msg << "C := " << alpha << " A' B + B' A" << beta << " C";
        Print( C, msg.str() );
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--r","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of update: N/T",'N');
        const Int m = Input("--m","height of result",100);
        const Int k = Input("--k","inner dimension",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );
        SetLocalTrr2kBlocksize<double>( nbLocal );
        SetLocalTrr2kBlocksize<Complex<double>>( nbLocal );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test Syr2k" << uploChar << transChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestSyr2k<double>( print, uplo, orientation, m, k, 3., 4., g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestSyr2k<Complex<double>>
        ( print, uplo, orientation, m, k, 
          Complex<double>(3), Complex<double>(4), g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
