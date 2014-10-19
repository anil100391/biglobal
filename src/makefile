CFLAGs	   = 
FFLAGS     = 
HOME       = /home/nebula/Dropbox/ME_project/biglobal/src
CPPFLAGS   = -I$(HOME)/include
FPPFLAGS   =
CLEANFILES = main
LIBD       = ../lib

include $(SLEPC_DIR)/conf/slepc_common

main: main.o chkopts
	${CLINKER} -o main main.o $(LIBD)/util.o -L$(LIBD)/ -lbiglobal ${SLEPC_EPS_LIB}
	rm *.o

all: cheb fd ndr eig
	ar q libbiglobal.a *.o
	ar q libbiglobal.a fd_grid.o

cheb:  
	$(PCC) -c cheb.c $(CCPPFLAGS) $(PCC_FLAGS)
	
fd:
	$(PCC) -c fd_grid/fd_grid.c $(CCPPFLAGS) $(PCC_FLAGS)
	$(PCC) -c fd.c $(CCPPFLAGS) $(PCC_FLAGS)

ndr:
	$(PCC) -c ndr.c $(CCPPFLAGS) $(PCC_FLAGS)
	
eig:
	$(PCC) -c eigen.c $(CCPPFLAGS) $(PCC_FLAGS)

install:
	mv libbiglobal.a ../lib/
	rm *.o
