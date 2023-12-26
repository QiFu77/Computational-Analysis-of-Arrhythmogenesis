all: SingleCell OneD TwoD

common = SingleCell/Cell.cc SingleCell/TP06.h SingleCell/TP06.cc

CC  	=	g++

CFLAGS	=	-w -O3 #-g:warning output to screen   -w:warning ignored

CFLAGS2	=	-fopenmp

CFLAGS3 =   -arch sm_61 -Xptxas -dlcm=cg

Initialization: $(common) Initialization.cc
	$(CC) $(CFLAGS) -o model_initialization $(common) Initialization.cc

SingleCell_TP06: $(common) SingleCell.cc
	$(CC) $(CFLAGS) -o model_SingleCell $(common) SingleCell.cc

SingleCell_ERP: $(common) SingleCell_ERP.cc
	$(CC) $(CFLAGS) -o model_SingleCell_ERP $(common) SingleCell_ERP.cc

OneD_VW: $(common) OneD_VW.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_oned_VW $(common) OneD_VW.cc

OneD: $(common) OneD.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_oned $(common) OneD.cc

TwoD: $(common) TwoD.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_twod $(common) TwoD.cc

ThreeD: $(common) ThreeD.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_threed $(common) ThreeD.cc

RestS1S2: $(common) RestCurveS1S2.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_rests1s2 $(common) RestCurveS1S2.cc

RestPCL: $(common) RestCurve_PCL.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_restpcl $(common) RestCurve_PCL.cc

clean:
	rm model_*
