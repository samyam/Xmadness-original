/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_ENTRYWISEMAP_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_PSEUDOSPECTRUM_INC

#include ELEM_FOXLI_INC
#include ELEM_DEMMEL_INC
#include ELEM_GRCAR_INC
#include ELEM_LOTKIN_INC
#include ELEM_UNIFORM_INC
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        Int r = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int matType =
            Input("--matType","0:uniform,1:Demmel,2:Lotkin,3:Grcar,4:FoxLi,"
                  "5:custom real,6:custom complex",1);
        bool quasi = Input("--quasi","Quasi-triang. real matrix?",true);
        const std::string basename =
            Input("--basename","basename of distributed Schur factor",
                  std::string("default"));
        const Int n = Input("--size","height of matrix",100);
        const Int nbAlg = Input("--nbAlg","algorithmic blocksize",96);
        const Real realCenter = Input("--realCenter","real center",0.);
        const Real imagCenter = Input("--imagCenter","imag center",0.);
        Real realWidth = Input("--realWidth","x width of image",0.);
        Real imagWidth = Input("--imagWidth","y width of image",0.);
        const Real numReal = Input("--numReal","num real chunks",2);
        const Real numImag = Input("--numImag","num imag chunks",2);
        const Int realSize = Input("--realSize","number of x samples",100);
        const Int imagSize = Input("--imagSize","number of y samples",100);
        const bool arnoldi = Input("--arnoldi","use Arnoldi?",true);
        const Int basisSize = Input("--basisSize","num Arnoldi vectors",10);
        const Int maxIts = Input("--maxIts","maximum pseudospec iter's",200);
        const Real psTol = Input("--psTol","tolerance for pseudospectra",1e-6);
        const Real uniformRealCenter = 
            Input("--uniformRealCenter","real center of uniform dist",0.);
        const Real uniformImagCenter =
            Input("--uniformImagCenter","imag center of uniform dist",0.);
        const Real uniformRadius =
            Input("--uniformRadius","radius of uniform dist",1.);
        const Int numBands = Input("--numBands","num bands for Grcar",3);
        const Real omega = Input("--omega","frequency for Fox-Li",16*M_PI);
        const bool progress = Input("--progress","print progress?",true);
        const bool deflate = Input("--deflate","deflate?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",false);
        const Int numSaveFreq = 
            Input("--numSaveFreq","numerical save frequency",-1);
        const Int imgSaveFreq = 
            Input("--imgSaveFreq","image save frequency",-1);
        const Int imgDispFreq =
            Input("--imgDispFreq","image display frequency",-1);
        const std::string numBase =
            Input("--numBase","numerical save basename",std::string("num"));
        const std::string imgBase =
            Input("--imgBase","image save basename",std::string("img"));
        const Int numFormatInt = Input("--numFormat","numerical format",2);
        const Int imgFormatInt = Input("--imgFormat","image format",8);
        const Int colorMapInt = Input("--colorMap","color map",0);
        const bool itCounts = Input("--itCounts","display iter. counts?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( mpi::Size(mpi::COMM_WORLD) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( mpi::COMM_WORLD, r, order );
        SetBlocksize( nbAlg );
        if( numFormatInt < 1 || numFormatInt >= FileFormat_MAX )
            LogicError("Invalid numerical format integer, should be in [1,",
                       FileFormat_MAX,")");
        if( imgFormatInt < 1 || imgFormatInt >= FileFormat_MAX )
            LogicError("Invalid image format integer, should be in [1,",
                       FileFormat_MAX,")");

        const FileFormat numFormat = static_cast<FileFormat>(numFormatInt);
        const FileFormat imgFormat = static_cast<FileFormat>(imgFormatInt);
        const ColorMap colorMap = static_cast<ColorMap>(colorMapInt);
        SetColorMap( colorMap );
        const C center(realCenter,imagCenter);
        const C uniformCenter(uniformRealCenter,uniformImagCenter);

        bool isReal = true;
        std::string matName;
        std::ostringstream os;
        DistMatrix<Real> AReal(g);
        DistMatrix<C> ACpx(g);
        switch( matType )
        {
        case 0: matName="uniform";
            Uniform( ACpx, n, n, uniformCenter, uniformRadius ); 
            MakeTriangular( UPPER, ACpx );
            isReal = false;
            break;
        case 1: matName="Demmel";
            Demmel( AReal, n ); 
            MakeTriangular( UPPER, AReal );
            isReal = true;
            break;
        case 2: matName="Lotkin";
            Lotkin( AReal, n );
            MakeTriangular( UPPER, AReal );
            isReal = true;
            break;
        case 3: matName="Grcar";
            Grcar( AReal, n, numBands ); 
            MakeTriangular( UPPER, AReal );
            isReal = true;
            break;
        case 4: matName="FoxLi";
            FoxLi( ACpx, n, omega );
            MakeTriangular( UPPER, ACpx );
            isReal = false;
            break;
        case 5: matName=basename;
            os << basename << "-" 
               << AReal.ColStride() << "x" << AReal.RowStride()
               << "-" << AReal.DistRank() << ".bin";
            AReal.Resize( n, n );
            read::Binary( AReal.Matrix(), os.str() );
            isReal = true;
            break;
        case 6: matName=basename;
            os << basename << "-" << ACpx.ColStride() << "x" << ACpx.RowStride()
               << "-" << ACpx.DistRank() << ".bin";
            ACpx.Resize( n, n );
            read::Binary( ACpx.Matrix(), os.str() );
            isReal = false;
            break;
        default:
            LogicError("Invalid matrix type");
        }
        if( display )
        {
            if( isReal )
                Display( AReal, "A" );
            else
                Display( ACpx, "A" );
        }
        if( write )
        {
            if( isReal )
            {
                Write( AReal, "A", numFormat );
                Write( AReal, "A", imgFormat );
            }
            else
            {
                Write( ACpx, "A", numFormat );
                Write( ACpx, "A", imgFormat );
            }
        }

        // Find a window if none is specified
        if( realWidth == 0. || imagWidth == 0. )
        {
            Real radius;
            if( isReal )
            {
                if( quasi )
                    radius = MaxNorm( schur::QuasiTriangEig( AReal ) );
                else
                    radius = MaxNorm( AReal.GetDiagonal() );
            }
            else
                radius = MaxNorm( ACpx.GetDiagonal() );
        
            const Real oneNorm = ( isReal ? OneNorm(AReal) : OneNorm(ACpx) );
            Real width;
            if( oneNorm == 0. && radius == 0. )
            {
                width = 1;
                if( mpi::WorldRank() == 0 )
                    std::cout << "Setting width to 1 to handle zero matrix"
                              << std::endl;
            }
            else if( radius >= 0.2*oneNorm )
            {
                width = 2.5*radius;
                if( mpi::WorldRank() == 0 )
                    std::cout << "Setting width to " << width
                              << " based on the spectral radius, " 
                              << radius << std::endl;
            }
            else
            {
                width = 0.8*oneNorm;
                if( mpi::WorldRank() == 0 )
                    std::cout << "Setting width to " << width
                              << " based on the one norm, " << oneNorm 
                              << std::endl;
            }
            realWidth = width;
            imagWidth = width;
        }

        PseudospecCtrl<Real> psCtrl;
        psCtrl.schur = true;
        psCtrl.maxIts = maxIts;
        psCtrl.tol = psTol;
        psCtrl.deflate = deflate;
        psCtrl.arnoldi = arnoldi;
        psCtrl.basisSize = basisSize;
        psCtrl.progress = progress;
        psCtrl.snapCtrl.imgSaveFreq = imgSaveFreq;
        psCtrl.snapCtrl.numSaveFreq = numSaveFreq;
        psCtrl.snapCtrl.imgDispFreq = imgDispFreq;
        psCtrl.snapCtrl.imgFormat = imgFormat;
        psCtrl.snapCtrl.numFormat = numFormat;
        psCtrl.snapCtrl.itCounts = itCounts;

        // Visualize/write the pseudospectrum within each window
        Timer timer;
        DistMatrix<Real> invNormMap(g);
        DistMatrix<Int> itCountMap(g);
        const Int xBlock = realSize / numReal;
        const Int yBlock = imagSize / numImag;
        const Int xLeftover = realSize - (numReal-1)*xBlock;
        const Int yLeftover = imagSize - (numImag-1)*yBlock;
        const Real xStep = realWidth/realSize;
        const Real yStep = imagWidth/imagSize;
        const C corner = center - C(realWidth/2,imagWidth/2);
        for( Int realChunk=0; realChunk<numReal; ++realChunk )
        {
            const Int realChunkSize = 
                ( realChunk==numReal-1 ? xLeftover : xBlock );
            const Real realChunkWidth = xStep*realChunkSize;
            for( Int imagChunk=0; imagChunk<numImag; ++imagChunk )
            {
                std::ostringstream chunkStream;
                chunkStream << "_" << realChunk << "_" << imagChunk;
                const std::string chunkTag = chunkStream.str();

                const Int imagChunkSize = 
                    ( imagChunk==numImag-1 ? yLeftover : yBlock );
                const Real imagChunkWidth = yStep*imagChunkSize;

                const C chunkCorner = corner + 
                    C(xStep*realChunk*xBlock,yStep*imagChunk*yBlock);
                const C chunkCenter = chunkCorner + 
                    0.5*C(xStep*realChunkSize,yStep*imagChunkSize);

                if( mpi::WorldRank() == 0 )
                    std::cout << "Starting computation for chunk centered at "
                              << chunkCenter << std::endl;
                mpi::Barrier( mpi::COMM_WORLD );
                timer.Start();
                psCtrl.snapCtrl.imgBase = matName+"-"+imgBase+chunkTag;
                psCtrl.snapCtrl.numBase = matName+"-"+numBase+chunkTag;
                if( isReal )
                {
                    if( quasi )
                        itCountMap = QuasiTriangularPseudospectrum
                        ( AReal, invNormMap, chunkCenter, 
                          realChunkWidth, imagChunkWidth, 
                          realChunkSize, imagChunkSize, psCtrl );
                    else
                        itCountMap = TriangularPseudospectrum
                        ( AReal, invNormMap, chunkCenter, 
                          realChunkWidth, imagChunkWidth, 
                          realChunkSize, imagChunkSize, psCtrl );
                }
                else
                {
                    itCountMap = TriangularPseudospectrum
                    ( ACpx, invNormMap, chunkCenter, 
                      realChunkWidth, imagChunkWidth, 
                      realChunkSize, imagChunkSize, psCtrl );
                }
                mpi::Barrier( mpi::COMM_WORLD );
                const double pseudoTime = timer.Stop();
                const Int numIts = MaxNorm( itCountMap );
                if( mpi::WorldRank() == 0 )
                {
                    std::cout << "num seconds=" << pseudoTime << "\n"
                              << "num iterations=" << numIts << std::endl;
                }
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
