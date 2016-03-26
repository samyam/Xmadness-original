/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/


/// \file examples/hello.cc
/// \brief Simplest example program for MADNESS
/// \defgroup hellowworldmad Hello world MADNESS style
/// \ingroup examples
///
/// Simplest program that initializes the MADNESS parallel runtime
/// using initialize(), makes a madness::World object, prints
/// a greeting, and then cleans up.
///
/// To initialize the MADNESS numerical environment you also need
/// \c startup(world,argc,argv) and should include mra/mra.h rather
/// than world/MADworld.h .

#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/FuseT/CopyOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/ReconstructOp.h>
#include <madness/mra/FuseT/NothingOp.h>
#include <madness/mra/FuseT/OpExecutor.h>

using namespace madness;;

static const double L		= 20;
static const long	k		= 8;
static const double thresh	= 1e-6;
//static const double thresh	= 1e-12;
static const double c		= 2.0;
static const double alpha	= 1.9; // Exponent

static double uinitial(const coord_3d& r) 
{
  const double x=r[0], y=r[1], z=r[2];
    return -2.0/(sqrt(x*x+y*y+z*z+1e-8));
    //return exp(-alpha*(2*x*x + y*y + z*z)) * pow(constants::pi/alpha, -1.5);
}

int main(int argc, char** argv) 
{
  initialize(argc,argv);
  World world(SafeMPI::COMM_WORLD);

  startup(world, argc, argv);

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_autorefine(false);
  FunctionDefaults<3>::set_cubic_cell(-L, L);

  if (world.rank() == 0)
  {
	  printf ("=================================================\n");
	  printf ("1. u0\n");
	  printf ("  1.1 compress by MADNESS\n");
	  printf ("  1.2 reconstruct by OpExecutor\n");
	  printf ("2. u1\n");
	  printf ("  2.1 compress by MADNESS\n");
	  printf ("  2.2 reconstruct by MADNESS\n");
	  printf ("=================================================\n");
  }

  real_function_3d u0 = real_factory_3d(world).f(uinitial);
  real_function_3d u1 = real_factory_3d(world).f(uinitial);
  u0.truncate();
  u1.truncate();

  double u1_norm	= u1.norm2();
  double u1_trace	= u1.trace();
  double u0_norm	= u0.norm2();
  double u0_trace	= u0.trace();

  if (world.rank() == 0) print("[u0] Initial norm", u0_norm,"trace", u0_trace);
  if (world.rank() == 0) print("[u1] Initial norm", u1_norm,"trace", u1_trace);
  world.gop.fence();


  // Make exponential of Vp
  real_function_3d result_factory = real_factory_3d(world);
  real_function_3d result(result_factory);
	
  if (world.rank() == 0)
  {
	  printf ("======================================\n");
	  printf ("Before u0 is compressed by MADNESS\n");
  }
  world.gop.fence();
  u0.compress();

  if (world.rank() == 0)
  {
	  printf ("======================================\n");
	  printf ("Before u0 is reconstructed by FuseT\n");
  }
  world.gop.fence();

  ReconstructOp<double,3> op1("Reconstruct",&result, &u0);
  OpExecutor<double,3> exe(world);
  exe.execute(&op1, false);
  world.gop.fence();

  if (world.rank() == 0)
  {
	  printf ("======================================\n");
	  printf ("After u0 is reconstructed by OpExecutor\n");
  } 
  world.gop.fence();
 
  double result_norm	= result.norm2();
  double result_trace	= result.trace();
  u0_norm	= u0.norm2();
  u0_trace	= u0.trace();


  if (world.rank() == 0) print("[u0] Result norm", u0_norm,"trace", u0_trace);
  world.gop.fence();
  if (world.rank() == 0) print("[result] Result norm", result_norm,"trace", result_trace);
  world.gop.fence();
/*
  //
  //	<---------
  //
  if (world.rank() == 0)
  {
	  printf ("======================================\n");
	  printf ("Before u1.compress()\n");
  }
  world.gop.fence();
  u1.compress();
  if (world.rank() == 0)
  {
	  printf ("======================================\n");
	  printf ("After u1.compress()\n");

	  printf ("======================================\n");
	  printf ("Before u1.reconstruct)\n");
  }
  world.gop.fence();
  u1.reconstruct();
  if (world.rank() == 0)
  {
	  printf ("======================================\n");
	  printf ("After u1.reconstructe()\n");
  }
  world.gop.fence();
  double result_norm = u1.norm2();
  double result_trace = u1.trace();
	
  if (world.rank() == 0) print("[u1] Result norm", result_norm," trace", result_trace);
  world.gop.fence();
*/	
  finalize();
  return 0;
}
