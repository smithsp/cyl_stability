FC=/usr/pppl/lff95/6.20c/bin/lf95 
FCOPTS=--dbl
LAPACK=/usr/pppl/lff95/6.20c/lapack-3.0
LOCAL=local
FEMOD=finite_elements_module
CYLFUNCS=cyl_funcs_module
VCYLMAT=vcyl_matrix_module
OBJ=${LOCAL}.o ${SPLINE}.o ${MODEL}.o model_funcs.o
LIBDIR=-L ${LAPACK} 
LIBS= -llapack -lblas 
NAG_LIBDIR = -L${NAG_ROOT}
NAG_LIBS = -lnag 
MODDIR = --mod ${NAG_ROOT}/nag_mod_dir
OUT=model.exe
TESTFE=test_finite_elements
OBJ=local.o sort_module.o ${FEMOD}.o ${CYLFUNCS}.o cyl.o
VOBJ=local.o sort_module.o ${FEMOD}.o ${VCYLMAT}.o vcyl.o
SRCDIR=src

local.o: ${SRCDIR}/${LOCAL}.f95
	${FC} $? ${FCOPTS} -c 
sort_module.o: ${SRCDIR}/sort_module.f95
	${FC} $? ${FCOPTS} -c 
${FEMOD}.o: ${SRCDIR}/${FEMOD}.f95
	${FC} $? ${FCOPTS} -c
${CYLFUNCS}.o: ${SRCDIR}/${CYLFUNCS}.f95
	${FC} $? ${FCOPTS} -c
${VCYLMAT}.o: ${SRCDIR}/${VCYLMAT}.f95
	${FC} $? ${FCOPTS} -c
${TESTFE}.o: ${SRCDIR}/${TESTFE}.f95
	${FC} ${FCOPTS} -c $?
vcyl.o: ${SRCDIR}/vcyl.f95
	${FC} $? ${FCOPTS} -c
cyl.o:  ${SRCDIR}/cyl.f95
	${FC} $? ${FCOPTS} -c
test_nag.o: ${SRCDIR}/test_nag.f95
	${FC} ${FCOPTS} -c $?
cyl: ${OBJ}
	${FC} ${OBJ} ${FCOPTS}  ${LIBDIR} ${LIBS} -o cyl.exe
vcyl: ${VOBJ}
	${FC} ${VOBJ} ${FCOPTS} -o vcyl.exe  ${LIBDIR} ${LIBS}  ${NAG_LIBDIR} ${NAG_LIBS}
testfe: local.o ${FEMOD}.o ${CYLFUNCS}.o ${TESTFE}.o
	${FC} local.o ${FEMOD}.o ${CYLFUNCS}.o ${TESTFE}.o ${FCOPTS} ${LIBDIR} ${LIBS}  -o testfe.exe
testnag: test_nag.o
	${FC} test_nag.o ${FCOPTS}  ${NAG_LIBDIR} ${NAG_LIBS}  -o testnag.exe 
clean:
	rm *.o *.mod *.exe
