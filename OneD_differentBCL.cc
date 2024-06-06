/*
 * General Code Structure (GCS) for one-dimensional simulation
 * 
 * A. Code sructure is re-organised using standard C++ to fit into my project.
 * B. The friendly style makes the project easier understood.
 * C. The version is more extendable whatever to further 1D,2D or another single cell model. 
 * 
 * Under Intellectual Property Protection.
 * 
 * 
 * Author      : Shugang Zhang <zhangshugang@hotmail.com>
 * Last update : 06-10-20023
 */

//***************************************注意如果用不同的单细胞模型需要注意几点：1.头文件的引用	2.后面对象的修改	3.输出文件的修改
//#include "SingleCell/TP06.h"
 #include "SingleCell/TPORd.cc"					//这一版本的ORd是TPORd模型


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "omp.h"

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENT
// #define HETEROGENEITY 

//***********************后期在做不同抗体的时候要修改下面的构造函数文件名称，以及输出的文件名称。  2、修改tp06.cc文件中的宏定义.  3.BCL 也要修改
//*****************************************说明：用来跑不同BCL下的ECG******
//*****************************************方法：不用读取已经initial的文件，采用new的方式，但是把numS1的数量增加到60次即让细胞达到稳态，同时只输出最后一个周期。
//*****************************************注意事项：1.每次需要修改TP06.cc的文件使它处在对应的细胞状态；     2.需要手动修改输出文件的文件名称。   3.修改对应的BCL；

//#define WT
//#define HETE
//#define HOMO
#define DOMINATE

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	// 1D parameters for all cell types
	double dx = 0.15; // mm
	double dt = 0.02; // ms
	int sanCellNum = 0;
	int atrialCellNum = 0;
	int epiCellNum = 0;
	int mCellNum = 0;
	int endoCellNum = 0;

	#if defined(VENT)	// 37 26 37    // 25 35 40
	endoCellNum = 25;
	mCellNum = 35;
	epiCellNum = 40;
	#endif
	// following parameters are for 1D or 2D.3D
	double atrialCoeff = 1.0*0.0195; // coefficient parameter in OneD strand, but needed to be validated by fitting CV.
	double sanCoeff = 0.1*0.0195; // coefficient parameter in OneD strand, but needed to be validated by fitting CV.
	double ventCoeff = 0.154; // 0.0106      0.154 for TP06          ****************?两个值的区别  在TPORd模型中用0.154

	// for ventricle
	#ifdef VENT
	int numS1 = 60;             //初始60		*************************达到稳态
	double BCL = 1000; // ms   325(2:1)   250(1:1)  1000(EAD) 500(normal) 750(normal)******************测试BCL 600-1500；*************************************************************************

	double stopTime = numS1*BCL; //ms
	double stimStrength = -52.0;//8.78; //8.78;//-8.78; // pA
	double stimDuration = 1.0;	// ms
	double stimStart = 0.0; // ms  // indicates the time point of beginning stimulus in a cycle
	#endif


	

	// --------start simulation--------	
	// CV calculation stuff
	double cvStartTime = 0;
	double cvEndTime = 0;
	int cvStartFlag = 0; // 0 for not start yet
	int cvEndFlag = 0;
	double cv = 0;

	// parallel stuff
	int coreNum = 8;//omp_get_num_procs();
	omp_set_num_threads(2 * coreNum);

	// strand initilization, diffusion stuff
	int cellNum = sanCellNum + atrialCellNum + epiCellNum + mCellNum + endoCellNum; // number of cells in OneD strand

	typedef Cell* CellPointer;
	//TP06* strand[cellNum]; // note that constructor contains the initializer
	TPORd* strand[cellNum];					//用的是TPORd模型

	double coeff[cellNum]; // diffusion parameters for each cell
	double dcoeff_dx[cellNum]; // first order derivative for each cell
	double oldV[cellNum];


	// assign coeff according to cell type
	for(int i = 0; i < cellNum; i++)
	{
		#ifdef VENT // set coeff to 'ventCoeff' whatever it was if VENT defined.
		// if(i == 59)
			// coeff[i] = 0.5*ventCoeff;
		if(i == 63)
			coeff[i] = 0.2*ventCoeff;
		// if(i == 61)
			// coeff[i] = 0.5*ventCoeff;
		else
			coeff[i] = ventCoeff;
		#endif
	}

	// Calculate the dcoeff/dx(i.e. dcoeff_dx in the code) in the 1D strand
	for(int i = 0; i < cellNum; i++)
	{
		if (i == 0) 
			dcoeff_dx[i] = (coeff[i+1] - coeff[i])/dx;
		else if (i == cellNum-1) 
			dcoeff_dx[i] = (coeff[i] - coeff[i-1])/dx;
		else
			dcoeff_dx[i] = (coeff[i+1] - coeff[i-1])/(2.0*dx);
	}


	#pragma omp parallel for schedule(static)

	for (int i = 0; i < cellNum; i++)
	{
		#ifdef WT
		FILE *initfile;
		if(i >= 0 && i < endoCellNum)
		{
			//strand[i] = new TP06(ENDO);
			strand[i] = new TPORd(ENDO);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
		// initfile = fopen("SingleCell/TP06InitialValues_WT_ENDO.dat","r");
		// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
		// fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			//strand[i] = new TP06(MCELL);
			strand[i] = new TPORd(MCELL);
			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
		// initfile = fopen("SingleCell/TP06InitialValues_WT_MCELL.dat","r");
		// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
		// fclose(initfile);
		}
		else // i < total cellnum
		{
			//strand[i] = new TP06(EPI);
			strand[i] = new TPORd(EPI);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
		// initfile = fopen("SingleCell/TP06InitialValues_WT_EPI.dat","r");
		// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
		// fclose(initfile);
		}	
		#endif
	#ifdef HETE
		FILE* initfile;
		if (i >= 0 && i < endoCellNum)
		{
			//strand[i] = new TP06(ENDO);
			strand[i] = new TPORd(ENDO);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_HETE_ENDO.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			//strand[i] = new TP06(MCELL);
			strand[i] = new TPORd(MCELL);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_HETE_MCELL.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
		else // i < total cellnum
		{
			//strand[i] = new TP06(EPI);
			strand[i] = new TPORd(EPI);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_HETE_EPI.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
	#endif

	#ifdef DOMINATE
		FILE* initfile;
		if (i >= 0 && i < endoCellNum)
		{
			//strand[i] = new TP06(ENDO);
			strand[i] = new TPORd(ENDO);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_ENDO.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			//strand[i] = new TP06(MCELL);
			strand[i] = new TPORd(MCELL);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_MCELL.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
		else // i < total cellnum
		{
			//strand[i] = new TP06(EPI);
			strand[i] = new TPORd(EPI);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_EPI.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
	#endif

	#ifdef HOMO
		FILE* initfile;
		if (i >= 0 && i < endoCellNum)
		{
			//strand[i] = new TP06(ENDO);
			strand[i] = new TPORd(ENDO);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_HOMO_ENDO.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
		else if (i < endoCellNum + mCellNum)
		{
			//strand[i] = new TP06(MCELL);
			strand[i] = new TPORd(MCELL);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_HOMO_MCELL.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
		else // i < total cellnum
		{
			//strand[i] = new TP06(EPI);
			strand[i] = new TPORd(EPI);

			// read in initial values (this is because the original init values is not stable yet)
			// if the initfile is not available, run the initialization.cc first
			// initfile = fopen("SingleCell/TP06InitialValues_HOMO_EPI.dat", "r");
			// strand[i]->readinAllStates(initfile);   // must be uncommented to get 0.719 m/s
			// fclose(initfile);
		}
	#endif



		// apply user configuration about dt
		strand[i]->setDt(dt);
	}
	//**********************结束初始化****************
	
	//需要手动修改文件名
	#ifdef WT
	//FILE *datafile = fopen("Outputs/VentOneDResultss_different-BCL_WT_1500.dat","w+");
	FILE *datafile = fopen("Outputs/VentOneDResultss_TPORd_different-BCL_WT_1500.dat","w+");
	#endif
	#ifdef HETE
		//FILE* datafile = fopen("Outputs/VentOneDResultss_different-BCL_HETE_.dat", "w+");
		FILE *datafile = fopen("Outputs/VentOneDResultss_TPORd_different-BCL_HETE_1500.dat","w+");
	#endif

	#ifdef DOMINATE
		//FILE* datafile = fopen("Outputs/VentOneDResults_different-BCL_DOMINATE_700.dat", "w+");
		//FILE *datafile = fopen("Outputs/VentOneDResultss_TPO6_numS1=60_different-BCL_DOMINATE_600.dat","w+");
		FILE *datafile = fopen("Outputs/VentOneDResultss_TPORd_different-BCL_DOMINATE_antibody60_1000.dat","w+");
	#endif

	#ifdef HOMO
		//FILE* datafile = fopen("Outputs/VentOneDResultss_different-BCL_HOMO_.dat", "w+");
		FILE *datafile = fopen("Outputs/VentOneDResultss_TPORd_different-BCL_HOMO_1000.dat","w+");
	#endif

       
	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time <= stopTime; time += dt, step++)
	{
		// 1. progress stats
		if(step%25000 == 0) // 25000 * dt ms = 0.125s 
			cout << "Progress = " << 100.0*time/stopTime << "\%." << endl;

		for(int i = 0; i < cellNum; i++)
		{
			oldV[i] = strand[i]->getV();
		}

		
		// 2. cell loop
		// #pragma omp parallel for schedule(static)
		for(int i = 0; i < cellNum; i++)
		{
			// ---stimulation applying or not, based on the time and cell location---
			// initialize all stim to 0.0
			strand[i]->setIstim(0.0);
			// Apply stimulus according to user configuration (about time)
			if(time - floor(time/BCL)*BCL >= stimStart && 
		   	   time - floor(time/BCL)*BCL < stimStart + stimDuration)
			{// first 5 cells get S1 stimulation
		    	if(i < 3 && i >= 0)
				{// cells get stimulation in certain duration
					// as strand[i] is a pointer, A->B should be used instead of A.B
					// ref. http://blog.csdn.net/qq457163027/article/details/54237782
					// cout << "here" << endl;
					strand[i]->setIstim(stimStrength);
				}
			}
			// following is S2 protocal, temporally commented out
			/**
			else if(time - (numS1-1)*BCL >= s2startTime && // s2 is an advanced stimulation, thus s2startTime ranges between 0~(BCL)
		   			time - (numS1-1)*BCL < s2startTime + stimDuration) 
			{// index from 10 to 15 cells get S2 stimulation
				if(i >= 10 && i < 16)
				// further constraint s2startTime to 0~(BCL - stimulation)
				{// cells get stimulation in certain duration
					// as strand[i] is a pointer, A->B should be used instead of A.B
					// ref. http://blog.csdn.net/qq457163027/article/details/54237782
					strand[i]->setIstim(s2stimStrength);
				}
			}
			**/


			// ---------calculate diffusion, i.e. dVgap---------
		
			double dVgap_dt = 0;
			double first_order;
			double second_order;

			// Step 1: calculate first and second order of membrane potential
			if(i == 0) 
			{
				// use strand[0] instead of "strand[-1]"
				first_order = (oldV[i+1] - oldV[i])/(1.0*dx);
				second_order = (oldV[i+1] + oldV[i] - 2.0*oldV[i])/(dx*dx);
			}
			else if(i > 0 && i < cellNum - 1) 
			{
				// normal case
				first_order = (oldV[i+1] - oldV[i-1])/(2.0*dx);
				second_order = (oldV[i+1] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}
			else if(i == cellNum - 1)
			{
				// use strand[cellNum-1] instead of "strand[cellNum]" as the latter is out of index
				first_order = (oldV[i] - oldV[i-1])/(1.0*dx);
				second_order = (oldV[i] + oldV[i-1] - 2.0*oldV[i])/(dx*dx);	
			}

			// Step 2: calculate dVgap according to equations
			dVgap_dt = dcoeff_dx[i]*first_order + coeff[i]*second_order;

			//if(step%1000 == 0) cout<<strand[i]->getVgap()<<endl;
			strand[i]->setDVgap_dt(dVgap_dt);
				
			// update
			strand[i]->update();
		}// end cell loop

		// 3. output file. Unfortunately, this part cannot run in paralell
	//	if( step % 10 == 0) //   输出所有周期
		if(floor(time / BCL) >= numS1 - 1 && step%10 == 0) // 50*dt = 1 ms once,output final cycle only 
		{
			for(int j = 0; j < cellNum; j++)
			{
				// write time before each time step
				if(j == 0)
					fprintf(datafile,"%4.10f\t", time);
				// write volt for each cell
				fprintf(datafile,"%4.10f\t", strand[j]->getV()); // unit: mV
				// write enter after each time step
				if(j == cellNum - 1)
					fprintf(datafile,"\n");
			}
		}// end Membrane Potential recording

		// 4. CV stats. Note that this can be only applied to non-selfpacing cells!!!
		// only record the last S1 data
		if (floor(time/BCL) == numS1-1)
		{ 
			// record propagation start and end time for CV calculation
			if(strand[10]->getV() >= -30 && cvStartFlag == 0)
			{
				cvStartTime = time;
				cout << "start = " << cvStartTime << endl;
				cvStartFlag = 1;
			}
			if(strand[90]->getV() >= -30 && cvEndFlag == 0)
			{
				cvEndTime = time;
				cout << "end = " << cvEndTime << endl;
				cvEndFlag = 1;
				cv = (dx * 80) / (cvEndTime - cvStartTime);
				cout << "duration = " << cvEndTime - cvStartTime << endl;
			}
		} // end CV calculation

	}// end of timeloop
	fclose(datafile);

	if(cvStartFlag == 1 && cvEndFlag == 1)
		cout << "CV = " << cv << " m/s." << endl;
	else
		cout << "Conduction failure!" << endl;
	printf("All done.\n");

	return 0;
}
