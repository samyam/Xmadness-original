/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_GEMM_INC
#include ELEM_TRSM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
using namespace std;
using namespace elem;

template<typename F> 
void TestTrsm
( bool print,
  LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  Int m, Int n, F alpha, const Grid& g )
{
    DistMatrix<F> A(g), X(g);

    if( side == LEFT )
        HermitianUniformSpectrum( A, m, 1, 10 );
    else
        HermitianUniformSpectrum( A, n, 1, 10 );
    auto S( A );
    MakeTriangular( uplo, S );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    Gemm( NORMAL, NORMAL, F(1)/alpha, S, X, Y );

    if( print )
    {
        Print( A, "A" );
        Print( S, "S" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Trsm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Trsm( side, uplo, orientation, diag, alpha, A, Y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 
        ( side==LEFT ? double(m)*double(m)*double(n)
                     : double(m)*double(n)*double(n) ) /(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. \n"
             << "  Time = " << runTime << " seconds. GFlops = " << gFlops 
             << endl;
    }
    if( print )
        Print( Y, "Y after solve" );
    Axpy( F(-1), X, Y );
    const auto SFrob = FrobeniusNorm( S );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    if( g.Rank() == 0 )
    {
        cout << "|| S ||_F = " << SFrob << "\n"
             << "|| X ||_F = " << XFrob << "\n"
             << "|| E ||_F = " << EFrob << "\n" << std::endl;
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
        const char sideChar = Input("--side","side to solve from: L/R",'L');
        const char uploChar = Input
            ("--uplo","lower or upper triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test Trsm" 
                 << sideChar << uploChar << transChar << diagChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestTrsm<double>( print, side, uplo, orientation, diag, m, n, 3., g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestTrsm<Complex<double>>
        ( print, side, uplo, orientation, diag, m, n, Complex<double>(3), g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
