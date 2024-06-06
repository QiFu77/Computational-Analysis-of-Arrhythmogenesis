/*
 * General Code Structure (GCS) for obtaining steady-state single cell initial values
 * 
 * A. Code sructure is re-organised using standard C++ to fit into my project.
 * B. The friendly style makes the project easier understood.
 * C. The version is more extendable whatever to further 1D,2D or another single cell model. 
 * 
 * Under Intellectual Property Protection.
 * 
 * 
 * Author      : Shugang Zhang <zhangshugang@hotmail.com>
 * Last update : 06-10-2023
 */

// #include "SingleCell/NeonatalRatAtria.cc"
// #include "SingleCell/RatSAN.cc"
// #include "SingleCell/GrandiCaMK.cc"
//#include "SingleCell/TP06.h"
// #include "SingleCell/ToRORD.cc"
 #include "SingleCell/TPORd.cc"  //****************************************我们的SingleCell中有两种单细胞模型，分别是TP06和TPORd，这里用的是TPORd。**********************

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENT

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.

	#ifdef VENT
	double BCL = 1000; // ms   //600-1500 
	double dt = 0.02; //ms
	double stopTime = 200*BCL; //ms
	double stimStrength = -52; // pA/pF
	double stimDuration = 1.0;	// ms
	double stimStart = 50.0; // ms  // indicates the time point of beginning stimulus in a cycle
	//typedef TP06 CellType;
	//typedef TP06* CellPointer;
	typedef TPORd CellType;
	typedef TPORd* CellPointer;
	#endif


	// --------start simulation--------
	//*************************************************************1.手动更改细胞类型	2.更改输出文件名字	3.输出文件要和TPORd.cc中细胞的类型相对应一定要检查。
	// note that constructor contains the initializer
	CellPointer cell = new CellType(MCELL);          //EPI   MCELL  ENDO
	#ifdef VENT
	//FILE *initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_EPI_antibody60.dat","w+");  //有抗体的文件名称规范
	FILE *initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_MCELL.dat","w+");  			 //没有抗体的文件名称规范
	#endif
	// apply user configuration about dt
	cell->setDt(dt);

	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%4000000 == 0) // 300000 * dt ms = 4.5s    0.005ms*4000000 = 20s
			cout << "Generating initial values, current progress = " << 100.0*time/stopTime << "\%." << endl;

		// Apply stimulus according to user configuration
		if(time - floor(time/BCL)*BCL >= stimStart && 
		   time - floor(time/BCL)*BCL < stimStart + stimDuration)	
		{
			cell->setIstim(stimStrength);
		}
		else
		{
			cell->setIstim(0.0);
		}

		// update all states and currents
		cell->update();

		// write file for each 1ms 
		if(time + dt > stopTime) // 50*dt = 1ms once
		{
			cell->outputAllStates(initfile);
			cout <<"Initial values of all states have been generated!" << endl;
		}
	}
	fclose(initfile);
	return 0;
}




