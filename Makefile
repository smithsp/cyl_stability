FC=/usr/pppl/lff95/6.20c/bin/lf95 
FCOPTS=--dbl 
LAPACK=/usr/pppl/lff95/6.20c/lapack-3.0
LOCAL=local
FEMOD=finite_elements_module
CYLFUNCS=cyl_funcs_module
VCYLFUNCS=vcyl_funcs_module
CYLMAT=cyl_matrix_module
VCYLMAT=vcyl_matrix_module
OBJ=${LOCAL}.o ${SPLINE}.o ${MODEL}.o model_funcs.o
LIBDIR=-L ${LAPACK} 
LIBS= -llapack -lblas 
NAG_LIBDIR = -L${NAG_ROOT}
NAG_LIBS = -lnag 
MODDIR = --mod ${NAG_ROOT}/nag_mod_dir
OUT=model.exe
TESTFE=test_finite_elements
OBJ=local.o sort_module.o ${FEMOD}.o ${CYLFUNCS}.o ${CYLMAT}.o cyl.o
VOBJ=local.o sort_module.o ${FEMOD}.o ${CYLFUNCS}.o ${VCYLFUNCS}.o ${CYLMAT}.o ${VCYLMAT}.o vcyl.o
EQOBJ=local.o sort_module.o ${FEMOD}.o ${VCYLMAT}.o 
SRCDIR=src

local.o: ${SRCDIR}/${LOCAL}.f95
	${FC} $? ${FCOPTS} -c 
sort_module.o: ${SRCDIR}/sort_module.f95
	${FC} $? ${FCOPTS} -c 
${FEMOD}.o: ${SRCDIR}/${FEMOD}.f95
	${FC} $? ${FCOPTS} -c
${CYLFUNCS}.o: ${SRCDIR}/${CYLFUNCS}.f95
	${FC} $? ${FCOPTS} -c
${VCYLFUNCS}.o: ${SRCDIR}/${VCYLFUNCS}.f95
	${FC} $? ${FCOPTS} -c
${CYLMAT}.o: ${SRCDIR}/${CYLMAT}.f95
	${FC} $? ${FCOPTS} -c
${VCYLMAT}.o: ${SRCDIR}/${VCYLMAT}.f95
	${FC} $? ${FCOPTS} -c ${NAG_LIBDIR} ${NAG_LIBS}
${TESTFE}.o: ${SRCDIR}/${TESTFE}.f95
	${FC} $? ${FCOPTS} -c
vcyl.o: ${SRCDIR}/vcyl.f95
	${FC} $? ${FCOPTS} -c
cyl.o:  ${SRCDIR}/cyl.f95
	${FC} $? ${FCOPTS} -c
test_nag.o: ${SRCDIR}/test_nag.f95
	${FC} $? ${FCOPTS} -c
plot_equilibrium.o: ${SRCDIR}/plot_equilibrium.f95
	${FC} $? ${FCOPTS} -c
test_eq.o: ${SRCDIR}/test_eq.f95
	${FC} $? ${FCOPTS} -c
cyl: ${OBJ}
	${FC} ${OBJ} ${FCOPTS}  ${LIBDIR} ${LIBS} -o cyl.exe ${NAG_LIBDIR} ${NAG_LIBS}
vcyl: ${VOBJ}
	${FC} ${VOBJ} ${FCOPTS} -o vcyl.exe  ${LIBDIR} ${LIBS}  ${NAG_LIBDIR} ${NAG_LIBS}
testfe: local.o ${FEMOD}.o ${CYLFUNCS}.o ${CYLMAT}.o ${TESTFE}.o sort_module.o
	${FC} local.o ${FEMOD}.o ${CYLFUNCS}.o ${CYLMAT}.o sort_module.o ${TESTFE}.o ${FCOPTS} ${LIBDIR} ${LIBS} ${NAG_LIBDIR} ${NAG_LIBS} -o testfe.exe
testnag: test_nag.o
	${FC} test_nag.o ${FCOPTS}  ${NAG_LIBDIR} ${NAG_LIBS}  -o testnag.exe 
plot_eq: ${EQOBJ} plot_equilibrium.o
	${FC}  ${EQOBJ} plot_equilibrium.o ${FCOPTS}  ${NAG_LIBDIR} ${NAG_LIBS}  -o plot_eq.exe 
test_eq: local.o ${CYLFUNCS}.o ${VCYLFUNCS}.o test_eq.o
	${FC} test_eq.o local.o ${CYLFUNCS}.o ${VCYLFUNCS}.o ${FCOPTS} -o test_eq.exe ${LIBDIR} ${LIBS}  ${NAG_LIBDIR} ${NAG_LIBS}
clean:
	rm *.o *.mod *.exe
