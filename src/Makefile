SHELL=/bin/bash
RELBASEDIR=..
RELBASEBLD=..
#----- Begin Boilerplate
ifneq (,$(wildcard ./Makefile))
# we are at the top level, in the source directory, phase 1 instance
# all build commands will be executed in an architecture specific subdirectory
.SUFFIXES:
# OBJDIR will be the build directiry for this architecture
ifeq (,$(EC_ARCH))
 OBJDIR := $(RELBASEBLD)/Build/base_arch/comp_arch$(SERIAL)
else
 OBJDIR := $(RELBASEBLD)/Build/$(EC_ARCH)$(SERIAL)
endif

MAKETARGET=$(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile \
	SRCDIR=$(CURDIR) $(MAKECMDGOALS)

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[[ -d $@ ]] || mkdir -p $@
	+@$(MAKETARGET)
	
Makefile : ;
%.mk :: ;
% :: $(OBJDIR) ; :

.PHONY: fullclean
fullclean:
	rm -rf $(EC_ARCH)

else
#----- End Boilerplate

# phase 2 instance is invoked in the build directory by phase 1

VPATH=$(SRCDIR)
# project building rules
include $(VPATH)/rpn_mpi.mk
#----- Begin Boilerplate
endif
#----- End Boilerplate
