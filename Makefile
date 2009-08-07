FC=/usr/pppl/lf64x/lf6480/bin/lfc
FCpar=/usr/pppl/lff95/8.x86_64-pkgs/mpich-1.2.7p1.x86_64/bin/mpif90
FCOPTS=--dbl #-g
#LAPACK=/usr/pppl/lff95/6.20c/lapack-3.0
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
VOBJ=local.o sort_module.o ${CYLFUNCS}.o ${FEMOD}.o  ${VCYLFUNCS}.o ${CYLMAT}.o ${VCYLMAT}.o vcyl.o
VPOBJ=local.o sort_module.o ${CYLFUNCS}.o ${FEMOD}.o  ${VCYLFUNCS}.o ${CYLMAT}.o ${VCYLMAT}.o vcyl_parallel.o
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
	${FC} $? ${FCOPTS} -c ${NAG_LIBDIR} ${NAG_LIBS}
vcyl_parallel.o: ${SRCDIR}/vcyl_parallel.f95
	${FCpar} $? ${FCOPTS} -c ${NAG_LIBDIR} ${NAG_LIBS}
cyl.o:  ${SRCDIR}/cyl.f95
	${FC} $? ${FCOPTS} -c
test_nag.o: ${SRCDIR}/test_nag.f95
	${FC} $? ${FCOPTS} -c ${MODDIR}
plot_equilibrium.o: ${SRCDIR}/plot_equilibrium.f95
	${FC} $? ${FCOPTS} -c
test_eq.o: ${SRCDIR}/test_eq.f95
	${FC} $? ${FCOPTS} -c
cyl: ${OBJ}
	${FC} ${OBJ} ${FCOPTS}  ${LIBDIR} ${LIBS} -o cyl.exe ${NAG_LIBDIR} ${NAG_LIBS}
vcyl: ${VOBJ}
	${FC} ${VOBJ} ${FCOPTS} -o vcyl.exe  ${LIBDIR} ${LIBS}  ${NAG_LIBDIR} ${NAG_LIBS}
vcylP: ${VPOBJ}
	${FCpar} ${VPOBJ} ${FCOPTS} -o vcyl_parallel.exe  ${LIBDIR} ${LIBS}  ${NAG_LIBDIR} ${NAG_LIBS}
testfe: local.o sort_module.o ${CYLFUNCS}.o ${FEMOD}.o ${CYLMAT}.o ${TESTFE}.o 
	${FC} local.o ${FEMOD}.o ${CYLFUNCS}.o sort_module.o ${CYLMAT}.o  ${TESTFE}.o ${FCOPTS} ${LIBDIR} ${LIBS} ${NAG_LIBDIR} ${NAG_LIBS} -o testfe.exe
testnag: local.o test_nag.o
	${FC} test_nag.o ${FCOPTS} ${MODDIR} ${NAG_LIBDIR} ${NAG_LIBS}  -o testnag.exe 
plot_eq: ${EQOBJ} plot_equilibrium.o
	${FC}  ${EQOBJ} plot_equilibrium.o ${FCOPTS}  ${NAG_LIBDIR} ${NAG_LIBS}  -o plot_eq.exe 
test_eq: local.o ${CYLFUNCS}.o ${VCYLFUNCS}.o test_eq.o
	${FC} test_eq.o local.o ${CYLFUNCS}.o ${VCYLFUNCS}.o ${FCOPTS} -o test_eq.exe ${LIBDIR} ${LIBS}  ${NAG_LIBDIR} ${NAG_LIBS}
clean:
	rm -f *.o *.mod *.exe
