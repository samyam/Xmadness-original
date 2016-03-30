#SUBDIRS = FusedTExamples
include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

# specify the git commit version in the info functions
GITREV="MADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\" "

# default location for basis sets etc
AM_CPPFLAGS += -D$(GITREV)


noinst_PROGRAMS = pno h2dft hedft hello hatom_energy h2 he tdse_example heat heat2 csqrt \
				  sdf_shape_tester test_gmres tdse1d vnucso nonlinschro sininteg functionio \
				  dataloadbal hatom_1d binaryop dielectric hehf 3dharmonic testsolver \
				  testspectralprop dielectric_external_field mp2 tiny oep h2dynamic newsolver \
				  tdhf cc2 nemo helium_exact density_smoothing siam_example gaussian \
<<<<<<< HEAD
				  copyEx compressEx reconstructEx innerEx fusedEx testFusedEx addEx multiplyEx derivativeEx

=======
				  copyEx compressEx reconstructEx innerEx fusedEx testFusedEx addEx multiplyEx \
				  matrixInnerEx
 
>>>>>>> 37e305ac5da2d05b850e1c264c4282e3b26bf609
if MADNESS_HAS_LIBXC
noinst_PROGRAMS += hefxc
endif

thisincludedir = $(includedir)/examples
thisinclude_HEADERS = molecularmask.h  nonlinsol.h spectralprop.h

AUTOMAKE_OPTIONS = subdir-objects

copyEx_SOURCES = FuseTExamples/CopyEx.cc
copyEx_LDADD = $(LIBWORLD) $(MRALIBS)

compressEx_SOURCES = FuseTExamples/compress.cc
compressEx_LDADD = $(LIBWORLD) $(MRALIBS)

<<<<<<< HEAD
derivativeEx_SOURCES = FuseTExamples/DerivativeEx.cc
derivativeEx_LDADD = $(LIBWORLD) $(MRALIBS)
=======
matrixInnerEx_SOURCES = FuseTExamples/MatrixInnerEx.cc
matrixInnerEx_LDADD = $(LIBWORLD) $(MRALIBS)
>>>>>>> 37e305ac5da2d05b850e1c264c4282e3b26bf609

reconstructEx_SOURCES = FuseTExamples/ReconstructEx.cc
reconstructEx_LDADD = $(LIBWORLD) $(MRALIBS)

innerEx_SOURCES = FuseTExamples/InnerEx.cc
innerEx_LDADD = $(LIBWORLD) $(MRALIBS)

multiplyEx_SOURCES = FuseTExamples/MultiplyEx.cc
multiplyEx_LDADD = $(LIBWORLD) $(MRALIBS)

addEx_SOURCES = FuseTExamples/AddEx.cc
addEx_LDADD = $(LIBWORLD) $(MRALIBS)

fusedEx_SOURCES = FuseTExamples/FusedEx.cc
fusedEx_LDADD = $(LIBWORLD) $(MRALIBS)

testFusedEx_SOURCES = FuseTExamples/Test_FusedEx.cc
testFusedEx_LDADD = $(LIBWORLD) $(MRALIBS)

pno_SOURCES = pno.cc 

siam_example_SOURCES = siam_example.cc 

mp2_SOURCES = mp2.cc 

helium_exact_SOURCES = helium_exact.cc 

oep_SOURCES = oep.cc 

tiny_SOURCES = tiny.cc

dataloadbal_SOURCES = dataloadbal.cc

dielectric_external_field_SOURCES = dielectric_external_field.cc

testspectralprop_SOURCES = testspectralprop.cc

testsolver_SOURCES = testsolver.cc

tdhf_SOURCES = tdhf.cc

cc2_SOURCES = cc2.cc

nemo_SOURCES = nemo.cc 

density_smoothing_SOURCES = density_smoothing.cc

newsolver_SOURCES = newsolver.cc

hedft_SOURCES = hedft.cc 
hedft_LDADD = $(LIBCHEM) $(MRALIBS)

h2dft_SOURCES = h2dft.cc
h2dft_LDADD = $(LIBCHEM) $(MRALIBS)

if MADNESS_HAS_LIBXC
hefxc_SOURCES = hefxc.cc 
hefxc_LDADD = $(LIBCHEM) $(MRALIBS)
endif

sininteg_SOURCES = sininteg.cc

3dharmonic_SOURCES = 3dharmonic.cc

functionio_SOURCES = functionio.cc

nonlinschro_SOURCES = nonlinschro.cc

heat_SOURCES = heat.cc

heat2_SOURCES = heat2.cc

hello_SOURCES = hello.cc

dielectric_SOURCES = dielectric.cc

hatom_energy_SOURCES = hatom_energy.cc

h2_SOURCES = h2.cc
 
h2dynamic_SOURCES = h2dynamic.cc

he_SOURCES = he.cc

hehf_SOURCES = hehf.cc

# h2hf_SOURCES = h2hf.cc

tdse_example_SOURCES = tdse_example.cc

tdse1d_SOURCES = tdse1d.cc

csqrt_SOURCES = csqrt.cc

sdf_shape_tester_SOURCES = sdf_shape_tester.cc

test_gmres_SOURCES = test_gmres.cc

vnucso_SOURCES = vnucso.cc

hatom_1d_SOURCES = hatom_1d.cc

binaryop_SOURCES = binaryop.cc

gaussian_SOURCES = gaussian.cc

LDADD = $(LIBCHEM) $(MRALIBS)

