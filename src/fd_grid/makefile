ALL: fd
CFLAGs	   = 
FFLAGS     = 
HOME       = /home/nebula/Dropbox/ME_project/biglobal/src
CPPFLAGS   = -I$(HOME)/include
FPPFLAGS   =
CLEANFILES = fd
#include $(PETSC_DIR)/conf/variables
#include $(PETSC_DIR)/conf/rules
include $(SLEPC_DIR)/conf/slepc_common

#fd: fd.o chkopts
#	mv fd.o ../lib/fd.o

fd:
	$(PCC) -c fd_grid.c $(CCPPFLAGS) $(PCC_FLAGS)
	mv fd_grid.o ../../lib/fd_grid.o

