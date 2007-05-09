FC=/usr/pppl/lff95/6.20c/bin/lf95 
FCOPTS= --dbl    
LAPACK=/usr/pppl/lff95/6.20c/lapack-3.0
LOCAL=local
FEMOD=finite_elements_module
CYLFUNCS=cyl_funcs_module
OBJ=${LOCAL}.o ${SPLINE}.o ${MODEL}.o model_funcs.o
LIBDIR=-L ${LAPACK}
LIBS= -llapack -lblas
OUT=model.exe
TESTFE=test_finite_elements
OBJFE=${CYLFUNCS}.o ${FEMOD}.o local.o 

local: ${LOCAL}.f95
	${FC} ${FCOPTS} -c $?
sort: sort_module.f95
	${FC} ${FCOPTS} -c $?
finite_elements: ${FEMOD}.f95
	${FC} ${FCOPTS} -c $?
cyl_funcs: ${CYLFUNCS}.f95
	${FC} ${FCOPTS} -c $?
cyl: cyl.f95
	make local
	make sort
	make finite_elements
	make cyl_funcs
	${FC} ${FCOPTS} cyl.f95 ${OBJFE} sort_module.o ${LIBDIR} ${LIBS} -o cyl.exe
	make clean
testfe: ${TESTFE}.f95
	make local
	make finite_elements
	make cyl_funcs
	${FC} ${FCOPTS}  $? ${OBJFE} ${LIBDIR} ${LIBS} -o testfe.exe
	make clean
clean:
	rm *.o *.mod
