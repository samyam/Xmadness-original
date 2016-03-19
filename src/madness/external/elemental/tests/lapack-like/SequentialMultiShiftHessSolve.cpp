/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKETRAPEZOIDAL_INC
#include ELEM_UPDATEDIAGONAL_INC
#include ELEM_AXPY_INC
#include ELEM_GEMM_INC
#include ELEM_MULTISHIFTHESSSOLVE_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_IDENTITY_INC
#include ELEM_UNIFORM_INC
using namespace elem;

// Test if (op(H) - mu_j I) x_j = y_j for each j.
// This is checked by testing the norm of  op(H) X - X Mu - Y.
template<typename F> 
void TestCorrectness
( UpperOrLower uplo, Orientation orientation, const Matrix<F>& H, 
  const Matrix<F>& shifts, const Matrix<F>& X, const Matrix<F>& Y,
  bool print, bool display )
{
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    
    Matrix<F> Z( Y );
    for( Int j=0; j<n; ++j )
    {
        auto x = LockedView( X, 0, j, m, 1 );
        auto z =       View( Z, 0, j, m, 1 );
        Axpy( shifts.Get(j,0), x, z );
    }

    Gemm( orientation, NORMAL, F(-1), H, X, F(1), Z );

    if( print && mpi::WorldRank() == 0 )
    {
        Print( H, "H" );
        Print( X, "X" );
        Print( Y, "Y" );
        Print( shifts, "shifts" );
        Print( Z, "-H X + X Mu + Y" );
    }
    if( display && mpi::WorldRank() == 0 )
    {
        Display( H, "H" );
        Display( X, "X" );
        Display( Y, "Y" );
        Display( shifts, "shifts" );
        Display( Z, "-H X + X Mu + Y" );
    }

    const Real YFrob = FrobeniusNorm( Y );
    const Real YInf = InfinityNorm( Y );
    const Real HFrob = FrobeniusNorm( H );
    const Real HInf = InfinityNorm( H );
    const Real ZFrob = FrobeniusNorm( Z );
    const Real ZInf = InfinityNorm( Z );
    if( mpi::WorldRank() == 0 )
    {
        std::cout << "    || H ||_F  = " << HFrob << "\n"
                  << "    || H ||_oo = " << HInf << "\n"
                  << "    || Y ||_F  = " << YFrob << "\n"
                  << "    || Y ||_oo = " << YInf << "\n"
                  << "    || H X - X Mu - Y ||_F  = " << ZFrob << "\n"
                  << "    || H X - X Mu - Y ||_oo = " << ZInf << "\n"
                  << std::endl;
    }
}

template<typename F>
void TestHessenberg
( UpperOrLower uplo, Orientation orientation, Int m, Int n, 
  bool testCorrectness, bool print, bool display )
{
    Matrix<F> H, X, Y, shifts;

    Uniform( H, m, m );
    UpdateDiagonal( H, F(5) ); // ensure that H-mu is far from zero
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );

    Uniform( X, m, n );
    Uniform( Y, m, n );
    Uniform( shifts, n, 1 );

    X = Y;
    if( mpi::WorldRank() == 0 )
    {
        std::cout << "  Starting Hessenberg solve...";
        std::cout.flush();
    }
    mpi::Barrier( mpi::COMM_WORLD );
    const double startTime = mpi::Time();
    MultiShiftHessSolve( uplo, orientation, F(1), H, shifts, X );
    mpi::Barrier( mpi::COMM_WORLD );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( mpi::WorldRank() == 0 )
    {
        std::cout << "DONE. " << std::endl
                  << "  Time = " << runTime << " seconds." << std::endl;
    }
    if( testCorrectness )
        TestCorrectness( uplo, orientation, H, shifts, X, Y, print, display );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    try
    {
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char orientChar = Input("--orient","orientation: N/T/C",'N');
        const Int m = Input("--m","height of Hessenberg matrix",100);
        const Int n = Input("--n","number of right-hand sides",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orient = CharToOrientation( orientChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            std::cout << "Double-precision:" << std::endl;
        TestHessenberg<double>
        ( uplo, orient, m, n, testCorrectness, print, display );

        if( commRank == 0 )
            std::cout << "Double-precision complex:" << std::endl;
        TestHessenberg<Complex<double>>
        ( uplo, orient, m, n, testCorrectness, print, display );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
