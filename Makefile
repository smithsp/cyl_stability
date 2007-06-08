FC=/usr/pppl/lff95/6.20c/bin/lf95 
FCOPTS= --dbl    
LAPACK=/usr/pppl/lff95/6.20c/lapack-3.0
LOCAL=local
FEMOD=finite_elements_module
CYLFUNCS=cyl_funcs_module
VCYLMAT=vcyl_matrix_module
OBJ=${LOCAL}.o ${SPLINE}.o ${MODEL}.o model_funcs.o
LIBDIR=-L ${LAPACK}
LIBS= -llapack -lblas
OUT=model.exe
TESTFE=test_finite_elements
OBJFE=${CYLFUNCS}.o ${FEMOD}.o local.o 
VOBJFE=${VCYLMAT}.o ${FEMOD}.o local.o 

local: src/${LOCAL}.f95
	${FC} ${FCOPTS} -c $?
sort: src/sort_module.f95
	${FC} ${FCOPTS} -c $?
finite_element: src/${FEMOD}.f95
	${FC} ${FCOPTS} -c $?
cyl_funcs: src/${CYLFUNCS}.f95
	${FC} ${FCOPTS} -c $?
vcyl_mat: src/${VCYLMAT}.f95
	${FC} ${FCOPTS} -c $?
cyl: src/cyl.f95
	make local
	make sort
	make finite_element
	make cyl_funcs
	${FC} ${FCOPTS} src/cyl.f95 ${OBJFE} sort_module.o ${LIBDIR} ${LIBS} -o cyl.exe
	make clean
vcyl: src/vcyl.f95
	make local
	make sort
	make finite_element
	make vcyl_mat
	${FC} ${FCOPTS} src/vcyl.f95 ${VOBJFE} sort_module.o ${LIBDIR} ${LIBS} -o vcyl.exe
	make clean
testfe: src/${TESTFE}.f95
	make local
	make finite_element
	make cyl_funcs
	${FC} ${FCOPTS}  $? ${OBJFE} ${LIBDIR} ${LIBS} -o testfe.exe
	make clean
clean:
	rm *.o *.mod
