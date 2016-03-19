/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/io.hpp"

#ifdef ELEM_HAVE_QT5

#include <QBoxLayout>

namespace elem {

SpyWindow::SpyWindow( QWidget* parent )
: QWidget(parent)
{
    DEBUG_ONLY(CallStackEntry cse("SpyWindow::SpyWindow"))
    matrix_ = 0;

    // For the real matrix
    QHBoxLayout* matrixLayout = new QHBoxLayout();
    spy_ = new SpyWidget();
    scroll_ = new QScrollArea();
    scroll_->setWidget( spy_ );
    matrixLayout->addWidget( scroll_ );
    setLayout( matrixLayout );

    setAttribute( Qt::WA_DeleteOnClose );

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

SpyWindow::~SpyWindow()
{ delete matrix_; }

void
SpyWindow::Spy( const Matrix<Int>* matrix, QString title )
{
    DEBUG_ONLY(CallStackEntry cse("SpyWindow::Spy"))
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    spy_->Spy( matrix );
}

} // namespace elem

#endif // ifdef ELEM_HAVE_QT5
