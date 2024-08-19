IGNORE:
CMPLR = gfortran
FFLAGS = -g -O0 -Wall -Wextra
LFLAGS = const_bse.h zdata.h 

METISSE_src_dir=$(BASE_dir)/METISSE/src
SSE_src_dir=$(BASE_dir)/SSE/src

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRC = \
test_bse.f hrdiag_remnant.f assign_remnant.f corerd.f comenv.f dgcore.f \
evolv2.f gntage.f instar.f kick.f mix.f mrenv.f ran3.f rl.f \
concatkstars.f comprad.f bpp_array.f checkstate.f deltat.f \
mlwind.f hrdiag.f star.f zcnsts.f assign_commons_COSMIC.f90

SSESRC = \
SSE_deltat.f \
SSE_mlwind.f SSE_hrdiag.f SSE_star.f \
SSE_zcnsts.f SSE_zfuncs.f SSE_gntage.f \


METISSESRC = \
METISSE_gntage.f90 METISSE_deltat.f90 \
METISSE_mlwind.f90 METISSE_hrdiag.f90 \
METISSE_star.f90 METISSE_zcnsts.f90 \
track_support.f90 z_support.f90 \
sse_support.f90  remnant_support.f90 \
interp_support.f90 comenv_lambda.f90 \
METISSE_miscellaneous.f90 

OBJT1 = $(SRC:.f=.o)

test:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o test
	mv test ../../debug/test

clean:
	rm -f *.o ../../debug/test