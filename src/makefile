COMPILER      = ${PCC}
FLAGS         = -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -march=native -mtune=native -std=c99

include ${PETSC_DIR}/lib/petsc/conf/variables

INC           = ${IGAP4_DIR}/include
LIB           = ${IGAP4_DIR}/lib
ODIR          = ${IGAP4_DIR}/lib/obj

_OBJ          = $(patsubst %.c,%.o,$(wildcard *.c))
OBJ           = $(patsubst %,$(ODIR)/%,$(_OBJ))

main: ${OBJ}

$(ODIR)/%.o: %.c
	${COMPILER} ${FLAGS} -I${INC} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -o $@ -c $<

