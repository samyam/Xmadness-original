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


  $Id$
*/

/// \file testbsh.cc
/// \brief test the bsh operator

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>

using namespace madness;

template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};


double aa;

double q(double r) {
    double val;
    if (r < 0.1e-4)
        val = (0.2e1 * exp(0.1e1 / aa / 0.4e1) * exp(-0.1e1 / aa / 0.4e1) * sqrt(aa) / sqrt(constants::pi) + exp(0.1e1 / aa / 0.4e1) * erf(0.1e1 / sqrt(aa) / 0.2e1) - exp(0.1e1 / aa / 0.4e1) + (0.2e1 / 0.3e1 * exp(0.1e1 / aa / 0.4e1) * exp(-0.1e1 / aa / 0.4e1) * (0.1e1 / 0.2e1 - aa) * sqrt(aa) / sqrt(constants::pi) + exp(0.1e1 / aa / 0.4e1) * erf(0.1e1 / sqrt(aa) / 0.2e1) / 0.6e1 - exp(0.1e1 / aa / 0.4e1) / 0.6e1) * r * r);
    else
        val = ((-exp((0.1e1 + 0.4e1 * aa * r) / aa / 0.4e1) + exp(-(-0.1e1 + 0.4e1 * aa * r) / aa / 0.4e1) + exp((0.1e1 + 0.4e1 * aa * r) / aa / 0.4e1) * erf((0.2e1 * aa * r + 0.1e1) / sqrt(aa) / 0.2e1) + exp(-(-0.1e1 + 0.4e1 * aa * r) / aa / 0.4e1) * erf((-0.1e1 + 0.2e1 * aa * r) / sqrt(aa) / 0.2e1)) / r / 0.2e1);

    return val / (4.0*constants::pi);
}


/// the result of the convolution
struct Qfunc : public FunctionFunctorInterface<double,3> {
    double operator()(const Vector<double,3>& x) const {
        double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        return q(r);
    }
};

template <typename T>
int test_bsh(World& world) {
    double mu = 1.0;
    std::vector<long> npt(3,201l);
    typedef Vector<double,3> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,3> > functorT;

    int success=0;
    if (world.rank() == 0)
        print("Test BSH operation, type =",
              archive::get_type_name<T>(),", ndim =",3);

    FunctionDefaults<3>::set_cubic_cell(-100,100);
    FunctionDefaults<3>::set_k(10);
    FunctionDefaults<3>::set_thresh(1e-5);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(true);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_truncate_on_project(false);

    const coordT origin(0.0);
    const double expnt = 100.0;
    aa = expnt;
    const double coeff = pow(expnt/constants::pi,1.5);

    // the input function to be convolved
    Function<T,3> f = FunctionFactory<T,3>(world).functor(functorT(new Gaussian<T,3>(origin, expnt, coeff)));
    f.truncate().reconstruct();

    double norm = f.trace();
    double ferr = f.err(Gaussian<T,3>(origin, expnt, coeff));
    if (world.rank() == 0) print("norm and error of the initial function", norm, ferr);


    // expnt=100 err=1e-9 use lo=2e-2 = .2/sqrt(expnt) and eps=5e-9

    // expnt=100 err=1e-7 use lo=2e-2 and eps=5e-7

    // expnt=100 err=1e-5 use lo=2e-e and eps=5e-5

    // expnt=100 err=1e-3 use lo=2e-2 and eps=5e-3


    SeparatedConvolution<T,3> op = BSHOperator<3>(world, mu, 1e-4, 1e-8);
    std::cout.precision(8);

    // apply the convolution operator on the input function f
    Function<T,3> ff = copy(f);
    if (world.rank() == 0) print("applying - 1");
    double start = cpu_time();
    Function<T,3> opf = op(ff);
    if (world.rank() == 0) print("done in time",cpu_time()-start);
    ff.clear();
    opf.verify_tree();
    double opferr = opf.err(Qfunc());
    if (world.rank() == 0) print("err in opf", opferr);
    if (world.rank() == 0) print("err in f", ferr);

    // here we are testing bsh, not the initial projection
    if ((opferr>ferr) and (opferr>FunctionDefaults<3>::get_thresh())) success++;

    return success;

    // FIXME: what comes here? Is it important??
    Function<double,3> qf = FunctionFactory<T,3>(world).functor(functorT(new Qfunc()));
    print("qf norm ", qf.norm2());
    print("opf norm", opf.norm2());

    opf.reconstruct();
//     for (int i=0; i<nn; ++i) {
//         double z=lo + i*range/double(nn-1);

//         double r = fabs(z)*sqrt(3.0);
//         coordT p(z);
// //         double r = z;
// //         coordT p(0.0); p[0] = z;

//         print(z, opf(p),q(r),opf(p)-q(r));
//     }

//     plotdx(opf, "opf.dx", FunctionDefaults<3>::get_cell(), npt);

    opf.truncate();
    Function<T,3> opinvopf = opf*(mu*mu);
    for (int axis=0; axis<3; ++axis) {
        //print("diffing",axis);
        //opinvopf.gaxpy(1.0,diff(diff(opf,axis),axis).compress(),-1.0);
    }

//     plotdx(opinvopf, "opinvopf.dx", FunctionDefaults<3>::get_cell(), npt);

    print("norm of (-del^2+mu^2)*G*f", opinvopf.norm2());
    Function<T,3> error = (f-opinvopf);
    print("error",error.norm2());

//     error.reconstruct();
//     plotdx(error, "err.dx", FunctionDefaults<3>::get_cell(), npt);


//     opinvopf.reconstruct();
//     f.reconstruct();
//     error.reconstruct();

// //     for (int i=0; i<101; ++i) {
// //         double z=-4 + 0.08*i;
// //         coordT p(z);
// //         print(z, opinvopf(p), f(p), opinvopf(p)/f(p), error(p));
// //     }

    opf.clear();
    opinvopf.clear();

    Function<T,3> g = (mu*mu)*f;
    for (int axis=0; axis<3; ++axis) {
        //g = g - diff(diff(f,axis),axis);
    }
    g = op(g);
    double derror=(g-f).norm2();
    print("norm of G*(-del^2+mu^2)*f",g.norm2());
    print("error",derror);
    if (derror> FunctionDefaults<3>::get_thresh()) success++;

    world.gop.fence();
    return success;

}


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);

    int success=0;
    try {
        startup(world,argc,argv);

        success=test_bsh<double>(world);

    }
    catch (const SafeMPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    finalize();

    return success;
}

