COMPILER      = ${PCC}
FLAGS         = -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -march=native -mtune=native -std=c99

include ${PETSC_DIR}/lib/petsc/conf/variables

INC           = ${IGAP4_DIR}/include
LIB           = ${IGAP4_DIR}/lib
ODIR          = ${IGAP4_DIR}/lib/obj

LIBNAME       = igap4

OBJ           = $(wildcard $(ODIR)/*.o)

main: 
	if [ ! -d ${IGAP4_DIR}/lib ]; then mkdir ${IGAP4_DIR}/lib; fi
	if [ ! -d ${IGAP4_DIR}/lib/obj ]; then mkdir ${IGAP4_DIR}/lib/obj; fi
	cd ${IGAP4_DIR}/src; make
	cd ${IGAP4_DIR}/mesh; make
	cd ${IGAP4_DIR}/phys; make
	cd ${IGAP4_DIR};

install:
	${COMPILER} ${FLAGS} -shared -Wl,-soname,lib${LIBNAME}.so.1 -o lib${LIBNAME}.so.1.0 ${OBJ} -lc
	mv lib${LIBNAME}.so.1.0 ${LIB}
	ln -sf ${LIB}/lib${LIBNAME}.so.1.0 ${LIB}/lib${LIBNAME}.so.1
	ln -sf ${LIB}/lib${LIBNAME}.so.1.0 ${LIB}/lib${LIBNAME}.so

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~




