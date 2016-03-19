/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_CONDITION_INC
#include ELEM_DETERMINANT_INC
#include ELEM_HILBERTSCHMIDT_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_NUCLEARNORM_INC
#include ELEM_TWONORM_INC
#include ELEM_HILBERT_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> H;
        Hilbert( H, n );
        if( display )
            Display( H, "Hilbert" );
        if( print )
            Print( H, "Hilbert matrix:" );

        // This is grossly inefficient due to recomputing the singular values
        // and Cholesky decomposition for several different operations, 
        // but it serves as an example of each function's usage
        const double cond = TwoCondition( H );
        const double det = HPDDeterminant( LOWER, H );
        const double hilbertSchmidt = HilbertSchmidt( H, H );
        const double twoNorm = HermitianTwoNorm( LOWER, H );
        const double frobNorm = HermitianFrobeniusNorm( LOWER, H );
        const double nuclearNorm = HermitianNuclearNorm( LOWER, H );

        if( mpi::WorldRank() == 0 )
        {
            std::cout << "kappa_2(H)   = " << cond << "\n"
                      << "det(H)       = " << det << "\n"
                      << "Tr(H' H)     = " << hilbertSchmidt << "\n"
                      << "|| H ||_F    = " << frobNorm << "\n"
                      << "|| H ||_*    = " << nuclearNorm << "\n"
                      << "|| H ||_2    = " << twoNorm << "\n"
                      << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
