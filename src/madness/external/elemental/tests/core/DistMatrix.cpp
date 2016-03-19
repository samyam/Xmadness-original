/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_IO_INC
#include ELEM_UNIFORM_INC
using namespace elem;

template<typename T,Dist AColDist,Dist ARowDist,Dist BColDist,Dist BRowDist>
void
Check( DistMatrix<T,AColDist,ARowDist>& A, 
       DistMatrix<T,BColDist,BRowDist>& B, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("Check"))
    const Grid& g = A.Grid();

    const Int commRank = g.Rank();
    const Int height = B.Height();
    const Int width = B.Width();
    DistMatrix<T,STAR,STAR> A_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> B_STAR_STAR(g);

    if( commRank == 0 )
    {
        std::cout << "Testing [" << DistToString(AColDist) << ","
                                 << DistToString(ARowDist) << "]"
                  << " <- ["     << DistToString(BColDist) << ","
                                 << DistToString(BRowDist) << "]...";
        std::cout.flush();
    }
    Int colAlign = SampleUniform<Int>(0,A.ColStride());
    Int rowAlign = SampleUniform<Int>(0,A.RowStride());
    mpi::Broadcast( colAlign, 0, mpi::COMM_WORLD );
    mpi::Broadcast( rowAlign, 0, mpi::COMM_WORLD );
    A.Align( colAlign, rowAlign );
    A = B;

    A_STAR_STAR = A;
    B_STAR_STAR = B;

    Int myErrorFlag = 0;
    for( Int j=0; j<width; ++j )
    {
        for( Int i=0; i<height; ++i )
        {
            if( A_STAR_STAR.GetLocal(i,j) != B_STAR_STAR.GetLocal(i,j) )
            {
                myErrorFlag = 1;
                break;
            }
        }
        if( myErrorFlag != 0 )
            break;
    }

    Int summedErrorFlag;
    mpi::AllReduce( &myErrorFlag, &summedErrorFlag, 1, mpi::SUM, g.Comm() );

    if( summedErrorFlag == 0 )
    {
        if( commRank == 0 )
            std::cout << "PASSED" << std::endl;
    }
    else
    {
        if( commRank == 0 )
            std::cout << "FAILED" << std::endl;
        if( print )
            Print( A, "A" );
        if( print ) 
            Print( B, "B" );
    }
}

template<typename T>
void
DistMatrixTest( Int m, Int n, const Grid& g, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixTest"))
    DistMatrix<T,MC,  MR  > A_MC_MR(g);
    DistMatrix<T,MC,  STAR> A_MC_STAR(g);
    DistMatrix<T,STAR,MR  > A_STAR_MR(g);
    DistMatrix<T,MR,  MC  > A_MR_MC(g);
    DistMatrix<T,MR,  STAR> A_MR_STAR(g);
    DistMatrix<T,STAR,MC  > A_STAR_MC(g);
    DistMatrix<T,VC,  STAR> A_VC_STAR(g);
    DistMatrix<T,STAR,VC  > A_STAR_VC(g);
    DistMatrix<T,VR,  STAR> A_VR_STAR(g);
    DistMatrix<T,STAR,VR  > A_STAR_VR(g);
    DistMatrix<T,STAR,STAR> A_STAR_STAR(g);

    // Communicate from A[MC,MR] 
    Uniform( A_MC_MR, m, n );
    Check( A_MC_STAR,   A_MC_MR, print );
    Check( A_STAR_MR,   A_MC_MR, print );
    Check( A_MR_MC,     A_MC_MR, print );
    Check( A_MR_STAR,   A_MC_MR, print );
    Check( A_STAR_MC,   A_MC_MR, print );
    Check( A_VC_STAR,   A_MC_MR, print );
    Check( A_STAR_VC,   A_MC_MR, print );
    Check( A_VR_STAR,   A_MC_MR, print );
    Check( A_STAR_VR,   A_MC_MR, print );
    Check( A_STAR_STAR, A_MC_MR, print );

    // Communicate from A[MC,*]
    Uniform( A_MC_STAR, m, n );
    Check( A_MC_MR,     A_MC_STAR, print );
    Check( A_STAR_MR,   A_MC_STAR, print );
    Check( A_MR_MC,     A_MC_STAR, print );
    Check( A_MR_STAR,   A_MC_STAR, print );
    Check( A_STAR_MC,   A_MC_STAR, print );
    Check( A_VC_STAR,   A_MC_STAR, print );
    Check( A_STAR_VC,   A_MC_STAR, print );
    Check( A_VR_STAR,   A_MC_STAR, print );
    Check( A_STAR_VR,   A_MC_STAR, print );
    Check( A_STAR_STAR, A_MC_STAR, print );

    // Communicate from A[*,MR]
    Uniform( A_STAR_MR, m, n );
    Check( A_MC_MR,     A_STAR_MR, print );
    Check( A_MC_STAR,   A_STAR_MR, print );
    Check( A_MR_MC,     A_STAR_MR, print );
    Check( A_MR_STAR,   A_STAR_MR, print );
    Check( A_STAR_MC,   A_STAR_MR, print );
    Check( A_VC_STAR,   A_STAR_MR, print );
    Check( A_STAR_VC,   A_STAR_MR, print );
    Check( A_VR_STAR,   A_STAR_MR, print );
    Check( A_STAR_VR,   A_STAR_MR, print );
    Check( A_STAR_STAR, A_STAR_MR, print );
    
    // Communicate from A[MR,MC]
    Uniform( A_MR_MC, m, n );
    Check( A_MC_MR,     A_MR_MC, print );
    Check( A_MC_STAR,   A_MR_MC, print );
    Check( A_STAR_MR,   A_MR_MC, print );
    Check( A_MR_STAR,   A_MR_MC, print );
    Check( A_STAR_MC,   A_MR_MC, print );
    Check( A_VC_STAR,   A_MR_MC, print );
    Check( A_STAR_VC,   A_MR_MC, print );
    Check( A_VR_STAR,   A_MR_MC, print );
    Check( A_STAR_VR,   A_MR_MC, print );
    Check( A_STAR_STAR, A_MR_MC, print );

    // Communicate from A[MR,*]
    Uniform( A_MR_STAR, m, n );
    Check( A_MC_MR,     A_MR_STAR, print );
    Check( A_MC_STAR,   A_MR_STAR, print );
    Check( A_STAR_MR,   A_MR_STAR, print );
    Check( A_MR_MC,     A_MR_STAR, print );
    Check( A_STAR_MC,   A_MR_STAR, print );
    Check( A_VC_STAR,   A_MR_STAR, print );
    Check( A_STAR_VC,   A_MR_STAR, print );
    Check( A_VR_STAR,   A_MR_STAR, print );
    Check( A_STAR_VR,   A_MR_STAR, print );
    Check( A_STAR_STAR, A_MR_STAR, print );

    // Communicate from A[*,MC]
    Uniform( A_STAR_MC, m, n );
    Check( A_MC_MR,     A_STAR_MC, print );
    Check( A_MC_STAR,   A_STAR_MC, print );
    Check( A_STAR_MR,   A_STAR_MC, print );
    Check( A_MR_MC,     A_STAR_MC, print );
    Check( A_MR_STAR,   A_STAR_MC, print );
    Check( A_VC_STAR,   A_STAR_MC, print );
    Check( A_STAR_VC,   A_STAR_MC, print );
    Check( A_VR_STAR,   A_STAR_MC, print );
    Check( A_STAR_VR,   A_STAR_MC, print );
    Check( A_STAR_STAR, A_STAR_MC, print );
 
    // Communicate from A[VC,*]
    Uniform( A_VC_STAR, m, n );
    Check( A_MC_MR,     A_VC_STAR, print );
    Check( A_MC_STAR,   A_VC_STAR, print );
    Check( A_STAR_MR,   A_VC_STAR, print );
    Check( A_MR_MC,     A_VC_STAR, print );
    Check( A_MR_STAR,   A_VC_STAR, print );
    Check( A_STAR_MC,   A_VC_STAR, print );
    Check( A_STAR_VC,   A_VC_STAR, print );
    Check( A_VR_STAR,   A_VC_STAR, print );
    Check( A_STAR_VR,   A_VC_STAR, print );
    Check( A_STAR_STAR, A_VC_STAR, print );

    // Communicate from A[*,VC]
    Uniform( A_STAR_VC, m, n );
    Check( A_MC_MR,     A_STAR_VC, print );
    Check( A_MC_STAR,   A_STAR_VC, print );
    Check( A_STAR_MR,   A_STAR_VC, print );
    Check( A_MR_MC,     A_STAR_VC, print );
    Check( A_MR_STAR,   A_STAR_VC, print );
    Check( A_STAR_MC,   A_STAR_VC, print );
    Check( A_VC_STAR,   A_STAR_VC, print );
    Check( A_VR_STAR,   A_STAR_VC, print );
    Check( A_STAR_VR,   A_STAR_VC, print );
    Check( A_STAR_STAR, A_STAR_VC, print );

    // Communicate from A[VR,*]
    Uniform( A_VR_STAR, m, n );
    Check( A_MC_MR,     A_VR_STAR, print );
    Check( A_MC_STAR,   A_VR_STAR, print );
    Check( A_STAR_MR,   A_VR_STAR, print );
    Check( A_MR_MC,     A_VR_STAR, print );
    Check( A_MR_STAR,   A_VR_STAR, print );
    Check( A_STAR_MC,   A_VR_STAR, print );
    Check( A_VC_STAR,   A_VR_STAR, print );
    Check( A_STAR_VC,   A_VR_STAR, print );
    Check( A_STAR_VR,   A_VR_STAR, print );
    Check( A_STAR_STAR, A_VR_STAR, print );

    // Communicate from A[*,VR]
    Uniform( A_STAR_VR, m, n );
    Check( A_MC_MR,     A_STAR_VR, print );
    Check( A_MC_STAR,   A_STAR_VR, print );
    Check( A_STAR_MR,   A_STAR_VR, print );
    Check( A_MR_MC,     A_STAR_VR, print );
    Check( A_MR_STAR,   A_STAR_VR, print );
    Check( A_STAR_MC,   A_STAR_VR, print );
    Check( A_VC_STAR,   A_STAR_VR, print );
    Check( A_STAR_VC,   A_STAR_VR, print );
    Check( A_VR_STAR,   A_STAR_VR, print );
    Check( A_STAR_STAR, A_STAR_VR, print );

    // Communicate from A[*,*]
    Uniform( A_STAR_STAR, m, n );
    Check( A_MC_MR,   A_STAR_STAR, print );
    Check( A_MC_STAR, A_STAR_STAR, print );
    Check( A_STAR_MR, A_STAR_STAR, print );
    Check( A_MR_MC,   A_STAR_STAR, print );
    Check( A_MR_STAR, A_STAR_STAR, print );
    Check( A_STAR_MC, A_STAR_STAR, print );
    Check( A_VC_STAR, A_STAR_STAR, print );
    Check( A_STAR_VC, A_STAR_STAR, print );
    Check( A_VR_STAR, A_STAR_STAR, print );
    Check( A_STAR_VR, A_STAR_STAR, print );
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
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print wrong matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );

        if( commRank == 0 )
            std::cout << "Testing with floats:" << std::endl;
        DistMatrixTest<float>( m, n, g, print );

        if( commRank == 0 )
            std::cout << "Testing with doubles:" << std::endl;
        DistMatrixTest<double>( m, n, g, print );

        if( commRank == 0 )
            std::cout << "Testing with single-precision complex:" << std::endl;
        DistMatrixTest<Complex<float>>( m, n, g, print );
        
        if( commRank == 0 )
            std::cout << "Testing with double-precision complex:" << std::endl;
        DistMatrixTest<Complex<double>>( m, n, g, print );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
