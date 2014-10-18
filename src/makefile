ALL: main
CFLAGs	   = 
FFLAGS     = 
HOME       = /home/nebula/Dropbox/ME_project/biglobal/src
CPPFLAGS   = -I$(HOME)/include
FPPFLAGS   =
CLEANFILES = main
LIBD       = ../lib

include $(SLEPC_DIR)/conf/slepc_common

main: main.o chkopts
	${CLINKER} -o main main.o fd.o $(LIBD)/fd_grid.o $(LIBD)/ndr.o $(LIBD)/eigen.o $(LIBD)/util.o ${SLEPC_EPS_LIB}

ndr:
	$(PCC) -c ndr.c $(CCPPFLAGS) $(PCC_FLAGS)
	mv ndr.o ../lib/ndr.o

eig:
	$(PCC) -c eigen.c $(CCPPFLAGS) $(PCC_FLAGS)
	mv eigen.o ../lib/eigen.o

fd:
	$(PCC) -c fd.c $(CCPPFLAGS) $(PCC_FLAGS)
