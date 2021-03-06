include $(VPATH)/RPN_MPI_version.inc
$(info this is RPN_MPI version $(RPN_MPI_version_s))

# general building rules
include $(VPATH)/Makefile.common

# sources specific (mechanically generated) dependencies
include $(VPATH)/dependencies.mk

LIB      = RPN_MPI$(SERIAL)

CLEAN    = RPN_MPI_fortran_stubs.F90 RPN_MPI_c_stubs.c mpi_stub.h \
           $(STUB_LIBRARY) $(LIBRARY) \
           $(VPATH)/rpn-comm_$(RPN_MPI_version_s)_multi.ssm \
           $(VPATH)/RPN_MPI_interfaces.hf $(VPATH)/dependencies.mk $(VPATH)/dependencies+.mk

CLEANDIRS= $(VPATH)/rpn-comm_$(RPN_MPI_version_s)_multi $(LIBDIR)

# build list of test programs and generate list of .Abs to build
TESTS    = $(shell ls $(VPATH)/TEST_*.[fF]90 | sed -e 's:.*/::' -e 's/[fF]90/Abs/g' | xargs)

FMODULES = RPN_MPI_mod.o
ifeq "$(SERIAL)" ""
MPI_VERSION = $(shell $(VPATH)/../tools/mpi_version.sh $(VTAG) )
else
MPI_VERSION = $(shell $(VPATH)/../tools/mpi_version.sh -none )
endif
LIBNAME  = $(LIB)_$(RPN_MPI_version)$(MPI_VERSION)
LIBRARY  = $(LIBDIR)/lib$(LIBNAME).a
STUB_LIBRARY = $(LIBDIR)/lib$(LIB)stubs_$(RPN_MPI_version)$(MPI_VERSION).a
SOURCES  = $(INCDECKS) $(CDECKS) $(FDECKS) $(HDECKS) $(F90DECKS)

DISTINCLUDES = $(VPATH)/RPN_MPI_interfaces.hf $(VPATH)/RPN_MPI.hf \
               $(VPATH)/RPN_MPI_types.inc $(VPATH)/RPN_MPI_ftoc.inc \
               $(VPATH)/RPN_MPI_mpi_layout.hf $(VPATH)/RPN_MPI_mpi_symbols.hf \
               $(VPATH)/RPN_MPI_mpi_definitions.hf $(VPATH)/RPN_MPI_mpif.hf \
               $(VPATH)/RPN_MPI_system_interfaces.hf \
               $(VPATH)/RPN_COMM_is_null.inc  $(VPATH)/RPN_COMM_constants.inc

ITF = $(VPATH)/RPN_MPI_interfaces.hf

lib: itf $(LIBRARY)

# OBJECTS comes from dependencies.mk
obj: $(OBJECTS)

all:   itf inc lib tests

full:  itf $(LIBRARY) $(TESTS) $(STUB_LIBRARY) $(LIBRARY).inc

tests:	$(TESTS)

stublib: $(STUB_LIBRARY)

# if producing a serial mpi library, add stubs to the library
ifeq "$(SERIAL)" ""
with-stubs: $(LIBRARY) stublib
	echo "WARNING: not in serial mode, stubs insertion will not be done"
else
with-stubs: $(LIBRARY)
	(cd $(VPATH)/../serial-mpi/stubs ; \
	rm -f *.o ; \
	s.cc -c -I../include RPN_MPI_c_stubs.c ; \
	s.f90 -c -I../include RPN_MPI_fortran_stubs.F90 ; \
	ar rcv $(LIBRARY) *.o ; \
	rm -f *.o ;)
endif

# special "library" that contains the include files
inc: $(LIBRARY).inc

itf: $(VPATH)/RPN_MPI_interfaces.hf

$(VPATH)/RPN_MPI_ptr.F90: $(VPATH)/../tools/gen_RPN_MPI_ptr.sh
	(cd $(VPATH) ; ../tools/gen_RPN_MPI_ptr.sh >$(VPATH)/RPN_MPI_ptr.F90)

#       extract InTfX minimalist interfaces first, then add InTf tagged interfaces
#       a .hf file is produced, included with #inclide <RPN_MPI_interfaces.hf>
#       #define IN_RPN_MPI_xxx 
#       allows to avoid conflict with functions/subroutines in source file RPN_MPI_xxx.F90
$(VPATH)/RPN_MPI_interfaces.hf: $(wildcard $(VPATH)/RPN_*.?90) $(wildcard $(VPATH)/RPN_*.c)
	(cd $(VPATH) ; rm -f RPN_MPI_interfaces.hf;)
	(cd $(VPATH) ; \
	for target in RPN_*.?90; \
	do grep -q '!InTfX!' $$target || continue ; \
	( echo "#if ! defined(IN_$${target%.*})" ; cat $$target | grep '!InTfX!' | sed 's/^!![ ]*/      /' ; echo "#endif" ) >>RPN_MPI_interfaces.hf ; \
	done;)
	(cd $(VPATH) ; \
	echo "interface" >> RPN_MPI_interfaces.hf ; \
	for target in RPN_*.?90 RPN_MPI_*.c ; \
	do ../tools/extract_interface.sh $$target >>RPN_MPI_interfaces.hf ; \
	done; \
	echo "end interface" >> RPN_MPI_interfaces.hf ; )

# (re)build dependencies using perl script rdedep.pl
$(VPATH)/dependencies.mk:
	rm -f $(VPATH)/dependencies.mk $(TMPDIR)/dependencies+.mk
	(cd $(VPATH) ; ../tools/rdedep.pl --flat_layout --out=$(TMPDIR)/dependencies+.mk $$(ls -1 RPN_* | grep -v TEST_)) || true
	mv $(TMPDIR)/dependencies+.mk $(VPATH)/dependencies.mk
	touch  $(VPATH)/dependencies.mk
#
dep_rm:
	rm -f $(VPATH)/dependencies.mk

dep: $(VPATH)/dependencies.mk

ssm-package:
	rm -rf $(VPATH)/rpn-comm_${RPN_MPI_version_s}_multi
	(cd $(VPATH) ; tar zxf ssmtemplate_1.0_all.ssm ; mv ssmtemplate_1.0_all rpn-comm_$(RPN_MPI_version_s)_multi )
	(cd $(VPATH) ; cp RPN_MPI_stubs.sh $(SOURCES) rpn-comm_$(RPN_MPI_version_s)_multi/src/.)
	(cd $(VPATH) ; cp Makefile Makefile.common Makefile.default Makefile.ECsetup *.mk \
	    rpn-comm_$(RPN_MPI_version_s)_multi/src/.)
	(cd $(VPATH) ; tar zcf rpn-comm_$(RPN_MPI_version_s)_multi.ssm rpn-comm_$(RPN_MPI_version_s)_multi)

RPN_MPI_fortran_stubs.o: $(VPATH)/RPN_MPI_stubs.sh
	$(SHELL) $(VPATH)/RPN_MPI_stubs.sh fortran

RPN_MPI_c_stubs.o: $(VPATH)/RPN_MPI_stubs.sh
	$(SHELL) $(VPATH)/RPN_MPI_stubs.sh c

$(STUB_LIBRARY): RPN_MPI_fortran_stubs.o RPN_MPI_c_stubs.o
	mkdir -p $(LIBDIR)
	ar rcv $(STUB_LIBRARY) RPN_MPI_fortran_stubs.o RPN_MPI_c_stubs.o
	(cd $(LIBDIR) ; ln -sf lib$(LIB)stubs_$(RPN_MPI_version)$(MPI_VERSION).a lib$(LIB)stubs$(MPI_VERSION).a)

# build static library and sorted list of objects in library, create link thst includes MPI_VERSION in name
$(LIBRARY): $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY)_ RPN_*.o
	mv $(LIBRARY)_ $(LIBRARY)
	ar t $(LIBRARY) | sort -u >$(VPATH)/objects.lst
	(cd $(LIBDIR) ; ln -sf lib$(LIBNAME).a  lib$(LIB)$(MPI_VERSION).a)

# check that contents of library (objects.lst) match reference list (REFERENCE.lst)
checkref:
	sort -u $(VPATH)/REFERENCE.lst >$(VPATH)/SORTED.lst
	diff $(VPATH)/SORTED.lst $(VPATH)/objects.lst

#$(VPATH)/includes: $(DISTINCLUDES)
$(LIBRARY).inc: $(DISTINCLUDES)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY).inc $(DISTINCLUDES)
	mkdir -p $(INCDIR)
	cp $(DISTINCLUDES) $(INCDIR)

tests_100 : itf $(LIBRARY)
	(cd $(VPATH) ; for i in TEST_1??.F90 ; do make $${i%.*}.Abs ; done ; )
