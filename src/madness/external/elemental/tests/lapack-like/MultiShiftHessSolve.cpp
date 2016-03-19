/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_CONJUGATE_INC
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
( UpperOrLower uplo, Orientation orientation, const DistMatrix<F,VC,STAR>& H, 
  const DistMatrix<F,VR,STAR>& shifts, 
  const DistMatrix<F,STAR,VR>& X, 
  const DistMatrix<F,STAR,VR>& Y, bool print, bool display )
{
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    auto modShifts( shifts );
    if( orientation == ADJOINT )
        Conjugate( modShifts );
    
    DistMatrix<F> Z( Y );
    for( Int j=0; j<n; ++j )
    {
        auto x = LockedView( X, 0, j, m, 1 );
        auto z =       View( Z, 0, j, m, 1 );
        Axpy( modShifts.Get(j,0), x, z );
    }
    {
        DistMatrix<F> H_MC_MR( H ), X_MC_MR(X);
        Gemm( orientation, NORMAL, F(-1), H_MC_MR, X_MC_MR, F(1), Z );
    }

    if( print )
    {
        Print( H, "H" );
        Print( X, "X" );
        Print( Y, "Y" );
        Print( shifts, "shifts" );
        Print( Z, "-H X + X Mu + Y" );
    }
    if( display )
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
  bool testCorrectness, bool print, bool display, const Grid& g )
{
    DistMatrix<F,VC,STAR> H(g);
    DistMatrix<F,STAR,VR> X(g), Y(g);
    DistMatrix<F,VR,STAR> shifts(g);

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
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
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

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid grid( comm, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orient = CharToOrientation( orientChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            std::cout << "Double-precision:" << std::endl;
        TestHessenberg<double>
        ( uplo, orient, m, n, testCorrectness, print, display, grid );

        if( commRank == 0 )
            std::cout << "Double-precision complex:" << std::endl;
        TestHessenberg<Complex<double>>
        ( uplo, orient, m, n, testCorrectness, print, display, grid );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
