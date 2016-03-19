/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_TRRK_LOCAL_HPP
#define ELEM_TRRK_LOCAL_HPP

#include ELEM_AXPYTRIANGLE_INC
#include ELEM_SCALETRAPEZOID_INC
#include ELEM_GEMM_INC

namespace elem {

namespace trrk {

#ifndef ELEM_RELEASE

void EnsureSame( const Grid& gA, const Grid& gB, const Grid& gC )
{
    if( gA != gB || gB != gC )
        LogicError("Grids must be the same");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,MC,STAR>& A, const DistMatrix<T>& C, std::string name )
{
    if( A.Height() != C.Height() || A.ColAlign() != C.ColAlign() )
        LogicError(name," not conformal with C");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,STAR,MC>& A, const DistMatrix<T>& C, std::string name )
{
    if( A.Width() != C.Height() || A.RowAlign() != C.ColAlign() )
        LogicError(name," not conformal with C");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,MR,STAR>& A, const DistMatrix<T>& C, std::string name )
{
    if( A.Height() != C.Width() || A.ColAlign() != C.RowAlign() )
        LogicError(name," not conformal with C");
}

template<typename T>
void EnsureConformal
( const DistMatrix<T,STAR,MR>& A, const DistMatrix<T>& C, std::string name )
{
    if( A.Width() != C.Width() || A.RowAlign() != C.RowAlign() )
        LogicError(name," not conformal with C");
}

template<typename T,Distribution UA,Distribution VA,
                    Distribution UB,Distribution VB>
void CheckInput
( const DistMatrix<T,UA,VA>& A, const DistMatrix<T,UB,VB>& B,
  const DistMatrix<T>& C )
{
    EnsureSame( A.Grid(), B.Grid(), C.Grid() );
    EnsureConformal( A, C, "A" );
    EnsureConformal( B, C, "B" );
}

// Local C := alpha A B + beta C
template<typename T>
void CheckInputNN( const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( A.Height() != C.Height() || B.Width()  != C.Width() ||
        A.Width()  != B.Height() || A.Height() != B.Width() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := alpha A B^{T/H} + beta C
template<typename T>
void CheckInputNT
( Orientation orientationOfB,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfB == NORMAL )
        LogicError("B must be (Conjugate)Transpose'd");
    if( A.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width()  != B.Width()  || A.Height() != B.Height() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := alpha A^{T/H} B + beta C
template<typename T>
void CheckInputTN
( Orientation orientationOfA,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfA == NORMAL )
        LogicError("A must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Width() != C.Width() ||
        A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

// Local C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void CheckInputTT
( Orientation orientationOfA,
  Orientation orientationOfB,
  const Matrix<T>& A, const Matrix<T>& B, const Matrix<T>& C )
{
    if( orientationOfA == NORMAL )
        LogicError("A must be (Conjugate)Transpose'd");
    if( orientationOfB == NORMAL )
        LogicError("B must be (Conjugate)Transpose'd");
    if( A.Width() != C.Height() || B.Height() != C.Width() ||
        A.Height() != B.Width() || A.Width() != B.Height() )
        LogicError
        ("Nonconformal LocalTrrk:\n",
         DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
}

#endif // ifndef ELEM_RELEASE

// Local C := alpha A B + beta C
template<typename T>
inline void
TrrkNNKernel
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("TrrkNNKernel");
        CheckInputNN( A, B, C );
    )
    Matrix<T> AT, AB;
    Matrix<T> BL, BR;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
    else
        Gemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );

    Gemm( NORMAL, NORMAL, alpha, AT, BL, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Gemm( NORMAL, NORMAL, alpha, AB, BR, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A B + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo, 
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrkKernel");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> AT(g), AB(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( NORMAL, NORMAL, alpha, AB, BL, T(1), CBL );
    else
        LocalGemm( NORMAL, NORMAL, alpha, AT, BR, T(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( NORMAL, NORMAL, alpha, AT, BL, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( NORMAL, NORMAL, alpha, AB, BR, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Local C := alpha A B^{T/H} + beta C
template<typename T>
inline void
TrrkNTKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("TrrkNTKernel");
        CheckInputNT( orientationOfB, A, B, C );
    )
    Matrix<T> AT, AB;
    Matrix<T> BT, BB;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
    else
        Gemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );

    Gemm( NORMAL, orientationOfB, alpha, AT, BT, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Gemm( NORMAL, orientationOfB, alpha, AB, BB, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A B^{T/H} + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrkKernel");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<T,MC,STAR> AT(g), AB(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionDown( A, AT, AB, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, T(1), CBL );
    else
        LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, T(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( NORMAL, orientationOfB, alpha, AT, BT, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( NORMAL, orientationOfB, alpha, AB, BB, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Local C := alpha A^{T/H} B + beta C
template<typename T>
inline void
TrrkTNKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("TrrkTNKernel");
        CheckInputTN( orientationOfA, A, B, C );
    )
    Matrix<T> AL, AR;
    Matrix<T> BL, BR;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
    else
        Gemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );

    Gemm( orientationOfA, NORMAL, alpha, AL, BL, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Gemm( orientationOfA, NORMAL, alpha, AR, BR, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrkKernel");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,STAR,MR> BL(g), BR(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionRight( B, BL, BR, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, T(1), CBL );
    else
        LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, T(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( orientationOfA, NORMAL, alpha, AL, BL, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( orientationOfA, NORMAL, alpha, AR, BR, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Local C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
inline void
TrrkTTKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("TrrkTTKernel");
        CheckInputTT( orientationOfA, orientationOfB, A, B, C );
    )
    Matrix<T> AL, AR;
    Matrix<T> BT, BB;
    Matrix<T> CTL, CTR,
              CBL, CBR;
    Matrix<T> DTL, DBR;

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        Gemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
    else
        Gemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );

    Gemm( orientationOfA, orientationOfB, alpha, AL, BT, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    Gemm( orientationOfA, orientationOfB, alpha, AR, BB, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
inline void
LocalTrrkKernel
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrkKernel");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    DistMatrix<T,STAR,MC> AL(g), AR(g);
    DistMatrix<T,MR,STAR> BT(g), BB(g);
    DistMatrix<T> CTL(g), CTR(g),
                  CBL(g), CBR(g);
    DistMatrix<T> DTL(g), DBR(g);

    const Int half = C.Height()/2;
    ScaleTrapezoid( beta, uplo, C );
    LockedPartitionRight( A, AL, AR, half );
    LockedPartitionDown( B, BT, BB, half );
    PartitionDownDiagonal
    ( C, CTL, CTR,
         CBL, CBR, half );

    if( uplo == LOWER )
        LocalGemm( orientationOfA, orientationOfB, alpha, AR, BT, T(1), CBL );
    else
        LocalGemm( orientationOfA, orientationOfB, alpha, AL, BB, T(1), CTR );

    DTL.AlignWith( CTL );
    LocalGemm( orientationOfA, orientationOfB, alpha, AL, BT, DTL );
    AxpyTriangle( uplo, T(1), DTL, CTL );

    DBR.AlignWith( CBR );
    LocalGemm( orientationOfA, orientationOfB, alpha, AR, BB, DBR );
    AxpyTriangle( uplo, T(1), DBR, CBR );
}

} // namespace trrk

namespace internal {

// Local C := alpha A B + beta C
template<typename T>
void TrrkNN
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrrkNN");
        CheckInputNN( A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkNNKernel( uplo, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AT, AB;
        Matrix<T> BL, BR;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( NORMAL, NORMAL, alpha, AB, BL, beta, CBL );
        else
            Gemm( NORMAL, NORMAL, alpha, AT, BR, beta, CTR );

        // Recurse
        TrrkNN( uplo, alpha, AT, BL, beta, CTL );
        TrrkNN( uplo, alpha, AB, BR, beta, CBR );
    }
}

// Local C := alpha A B^{T/H} + beta C
template<typename T>
void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrrkNT");
        CheckInputNT( orientationOfB, A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkNTKernel( uplo, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AT, AB;
        Matrix<T> BT, BB;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( NORMAL, orientationOfB, alpha, AB, BT, beta, CBL );
        else
            Gemm( NORMAL, orientationOfB, alpha, AT, BB, beta, CTR );

        // Recurse
        TrrkNT( uplo, orientationOfB, alpha, AT, BT, beta, CTL );
        TrrkNT( uplo, orientationOfB, alpha, AB, BB, beta, CBR );
    }
}

// Local C := alpha A^{T/H} B + beta C
template<typename T>
void TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrrkTN");
        CheckInputTN( orientationOfA, A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkTNKernel( uplo, orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AL, AR;
        Matrix<T> BL, BR;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( orientationOfA, NORMAL, alpha, AR, BL, beta, CBL );
        else
            Gemm( orientationOfA, NORMAL, alpha, AL, BR, beta, CTR );

        // Recurse
        TrrkTN( uplo, orientationOfA, alpha, AL, BL, beta, CTL );
        TrrkTN( uplo, orientationOfA, alpha, AR, BR, beta, CBR );
    }
}

// Local C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrrkTT");
        CheckInputTT( orientationOfA, orientationOfB, A, B, C );
    )
    if( C.Height() < LocalTrrkBlocksize<T>() )
    {
        TrrkTTKernel
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        Matrix<T> AL, AR;
        Matrix<T> BT, BB;
        Matrix<T> CTL, CTR,
                  CBL, CBR;

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            Gemm( orientationOfA, orientationOfB, alpha, AR, BT, beta, CBL );
        else
            Gemm( orientationOfA, orientationOfB, alpha, AL, BB, beta, CTR );

        // Recurse
        TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, AL, BT, beta, CTL );
        TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, AR, BB, beta, CBR );
    }
}

} // namespace internal

// Distributed C := alpha A B + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,  STAR>& A,
           const DistMatrix<T,STAR,MR  >& B,
  T beta,        DistMatrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrk");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,MC,STAR> AT(g), AB(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( NORMAL, NORMAL, alpha, AB, BL, beta, CBL );
        else
            LocalGemm( NORMAL, NORMAL, alpha, AT, BR, beta, CTR );

        // Recurse
        LocalTrrk( uplo, alpha, AT, BL, beta, CTL );
        LocalTrrk( uplo, alpha, AB, BR, beta, CBR );
    }
}

// Distributed C := alpha A B^{T/H} + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A,
           const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrk");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,MC,STAR> AT(g), AB(g);
        DistMatrix<T,MR,STAR> BT(g), BB(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionDown( A, AT, AB, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( NORMAL, orientationOfB, alpha, AB, BT, beta, CBL );
        else
            LocalGemm( NORMAL, orientationOfB, alpha, AT, BB, beta, CTR );

        // Recurse
        LocalTrrk( uplo, orientationOfB, alpha, AT, BT, beta, CTL );
        LocalTrrk( uplo, orientationOfB, alpha, AB, BB, beta, CBR );
    }
}

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A,
           const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrk");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel( uplo, orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,STAR,MR> BL(g), BR(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionRight( B, BL, BR, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm( orientationOfA, NORMAL, alpha, AR, BL, beta, CBL );
        else
            LocalGemm( orientationOfA, NORMAL, alpha, AL, BR, beta, CTR );

        // Recurse
        LocalTrrk( uplo, orientationOfA, alpha, AL, BL, beta, CTL );
        LocalTrrk( uplo, orientationOfA, alpha, AR, BR, beta, CBR );
    }
}

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC  >& A,
           const DistMatrix<T,MR,  STAR>& B,
  T beta,        DistMatrix<T>& C )
{
    using namespace trrk;
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrrk");
        CheckInput( A, B, C );
    )
    const Grid& g = C.Grid();

    if( C.Height() < g.Width()*LocalTrrkBlocksize<T>() )
    {
        LocalTrrkKernel
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on corner
        // and recurse on CTL and CBR.
        DistMatrix<T,STAR,MC> AL(g), AR(g);
        DistMatrix<T,MR,STAR> BT(g), BB(g);
        DistMatrix<T> CTL(g), CTR(g),
                      CBL(g), CBR(g);

        const Int half = C.Height() / 2;
        LockedPartitionRight( A, AL, AR, half );
        LockedPartitionDown( B, BT, BB, half );
        PartitionDownDiagonal
        ( C, CTL, CTR,
             CBL, CBR, half );

        if( uplo == LOWER )
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AR, BT, beta, CBL );
        else
            LocalGemm
            ( orientationOfA, orientationOfB, alpha, AL, BB, beta, CTR );

        // Recurse
        LocalTrrk
        ( uplo, orientationOfA, orientationOfB, alpha, AL, BT, beta, CTL );
        LocalTrrk
        ( uplo, orientationOfA, orientationOfB, alpha, AR, BB, beta, CBR );
    }
}

} // namespace elem

#endif // ifndef ELEM_TRRK_LOCAL_HPP
