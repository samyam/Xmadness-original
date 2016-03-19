/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HEMM_INC
#include ELEM_TRMM_INC
#include ELEM_TWOSIDEDTRMM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( bool print, UpperOrLower uplo, UnitOrNonUnit diag,
  const DistMatrix<F>& A, const DistMatrix<F>& B, const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    const Int k=100;
    DistMatrix<F> X(g), Y(g), Z(g);
    Uniform( X, m, k );
    Y = X;
    Zeros( Z, m, k );

    if( uplo == LOWER )
    {
        // Test correctness by comparing the application of A against a 
        // random set of k vectors to the application of 
        // tril(B)^H AOrig tril(B)
        Trmm( LEFT, LOWER, NORMAL, diag, F(1), B, Y );
        Hemm( LEFT, LOWER, F(1), AOrig, Y, F(0), Z );
        Trmm( LEFT, LOWER, ADJOINT, diag, F(1), B, Z );
        Hemm( LEFT, LOWER, F(-1), A, X, F(1), Z );
        Real infNormOfAOrig = HermitianInfinityNorm( uplo, AOrig );
        Real frobNormOfAOrig = HermitianFrobeniusNorm( uplo, AOrig );
        Real infNormOfA = HermitianInfinityNorm( uplo, A );
        Real frobNormOfA = HermitianFrobeniusNorm( uplo, A );
        Real oneNormOfError = OneNorm( Z );
        Real infNormOfError = InfinityNorm( Z );
        Real frobNormOfError = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
        {
            cout << "||AOrig||_1 = ||AOrig||_oo     = "
                 << infNormOfAOrig << "\n"
                 << "||AOrig||_F                    = "
                 << frobNormOfAOrig << "\n"
                 << "||A||_1 = ||A||_oo             = "
                 << infNormOfA << "\n"
                 << "||A||_F                        = "
                 << frobNormOfA << "\n"
                 << "||A X - L^H AOrig L X||_1  = "
                 << oneNormOfError << "\n"
                 << "||A X - L^H AOrig L X||_oo = " 
                 << infNormOfError << "\n"
                 << "||A X - L^H AOrig L X||_F  = "
                 << frobNormOfError << endl;
        }
    }
    else
    {
        // Test correctness by comparing the application of A against a 
        // random set of k vectors to the application of 
        // triu(B) AOrig triu(B)^H
        Trmm( LEFT, UPPER, ADJOINT, diag, F(1), B, Y );
        Hemm( LEFT, UPPER, F(1), AOrig, Y, F(0), Z );
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), B, Z );
        Hemm( LEFT, UPPER, F(-1), A, X, F(1), Z );
        Real infNormOfAOrig = HermitianInfinityNorm( uplo, AOrig );
        Real frobNormOfAOrig = HermitianFrobeniusNorm( uplo, AOrig );
        Real infNormOfA = HermitianInfinityNorm( uplo, A );
        Real frobNormOfA = HermitianFrobeniusNorm( uplo, A );
        Real oneNormOfError = OneNorm( Z );
        Real infNormOfError = InfinityNorm( Z );
        Real frobNormOfError = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
        {
            cout << "||AOrig||_1 = ||AOrig||_oo     = "
                 << infNormOfAOrig << "\n"
                 << "||AOrig||_F                    = "
                 << frobNormOfAOrig << "\n"
                 << "||A||_1 = ||A||_oo             = "
                 << infNormOfA << "\n"
                 << "||A||_F                        = "
                 << frobNormOfA << "\n"
                 << "||A X - U AOrig U^H X||_1  = "
                 << oneNormOfError << "\n"
                 << "||A X - U AOrig U^H X||_oo = " 
                 << infNormOfError << "\n"
                 << "||A X - U AOrig U^H X||_F  = "
                 << frobNormOfError << endl;
        }
    }
}

template<typename F> 
void TestTwoSidedTrmm
( bool testCorrectness, bool print, UpperOrLower uplo, UnitOrNonUnit diag, 
  Int m, const Grid& g )
{
    DistMatrix<F> A(g), B(g), AOrig(g);

    Zeros( A, m, m );
    Zeros( B, m, m );
    MakeHermitianUniformSpectrum( A, 1, 10 );
    MakeHermitianUniformSpectrum( B, 1, 10 );
    MakeTriangular( uplo, B );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting reduction to Hermitian standard EVP...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TwoSidedTrmm( uplo, diag, A, B );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    double gFlops = Pow(double(m),3.)/(runTime*1.e9);
    if( IsComplex<F>::val )
        gFlops *= 4.;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = "
             << gFlops << endl;
    }
    if( print )
        Print( A, "A after reduction" );
    if( testCorrectness )
        TestCorrectness( print, uplo, diag, A, B, AOrig );
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
        const char uploChar = Input
            ("--uplo","lower or upper triangular storage: L/U",'L');
        const char diagChar = Input("--unit","(non-)unit diagonal: N/U",'N');
        const Int m = Input("--m","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test TwoSidedTrmm" << uploChar << diagChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestTwoSidedTrmm<double>( testCorrectness, print, uplo, diag, m, g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestTwoSidedTrmm<Complex<double>>
        ( testCorrectness, print, uplo, diag, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
