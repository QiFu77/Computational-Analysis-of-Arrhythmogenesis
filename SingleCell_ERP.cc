/*
 * General Code Structure (GCS) for three-dimensional simulation
 *
 * A. Code sructure is re-organised using standard C++ to fit into my project.
 * B. The friendly style makes the project easier understood.
 * C. The version is more extendable whatever to further 1D,2D or another single cell model.
 *
 * Under Intellectual Property Protection.
 *
 *
 * Author      : Shugang Zhang <zhangshugang@hotmail.com>
 * Last update : 11-12-2023
 */


//#include "SingleCell/TP06.h"
 #include "SingleCell/TPORd.cc"

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENTRICLE
#define EPSILON 1e-7
//**************1.选择构造文件细胞，注意有两处，S1和S2两处地方都要选择********2.下面还需要修改输出文件的名字。3.修改tp06.cc里面的文件，确认是处于病理状态还是药理状态。
#define ENDO1
 //#define EPI1
//#define MCELL1

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	#if defined(ATRIA) || defined(VENTRICLE)
	double BCL = 1000; // ms   // NOTE THERE ARE TWO INIT FILES TO BE CHANGED!!!
	double numS1 = 2; // only 2 is needed cause cell have been initialized using pre-trained file (initialization.dat)
	double dt = 0.02; //ms
	double stopTime = numS1*BCL + BCL; //ms, another BCL is added for calculating S2 and ERP.
	double stimStrength = -52; // pA/pF   -6.0pA/pF(-0.6nA) for RatAtrial // -12.5pA/pF for CaMKII; -8.78pA for neonatalRatAtrial
	double stimDuration = 1;   // ms
	double stimStart = 50.0;   // ms  // indicates the time point of beginning stimulus in a cycle
	// typedef GrandiCaMKII CellType;
	// typedef GrandiCaMKII* CellPointer;
	// typedef ORdHumanVentricle CellType;
	// typedef ORdHumanVentricle* CellPointer;
	// typedef NeonatalRatAtria CellType;
	// typedef NeonatalRatAtria* CellPointer;
	
	//****************************************我们的SingleCell中有两种单细胞模型，分别是TP06和TPORd，这里用的是TPORd。**********************
	//typedef TP06 CellType;
	//typedef TP06* CellPointer;
	typedef TPORd CellType;
	typedef TPORd* CellPointer;

	// statistics for single cell
	double apd20;
	double apd25;
	double apd50;
	double apd80;
	double apd90;
	double dvdt_max;
	double rest_potential;
	double amplitude;
	double overshoot;
	#endif




	// --------start simulation--------

	// initialize from file
	FILE *initfile;
	

	#ifdef EPI1
	CellPointer cell = new CellType(EPI);
	//initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_EPI_antibody60.dat","r"); // NOTE THERE ARE TWO INIT FILES TO BE CHANGED!!!
	initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_EPI.dat","r");
	cell->readinAllStates(initfile);
	#endif 

	#ifdef ENDO1
	CellPointer cell = new CellType(ENDO);
	//initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_ENDO_antibody60.dat","r");
	initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_ENDO.dat","r");
	cell->readinAllStates(initfile);
	#endif 

	#ifdef MCELL1
	CellPointer cell = new CellType(MCELL);
	//initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_MCELL_antibody60.dat","r");
	initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_MCELL.dat","r");
	cell->readinAllStates(initfile);
	#endif 

	#ifdef VENTRICLE
	//FILE *datafile = fopen("Outputs/VentERP_TPORd_DOMINATE_ENDO.dat","w+");
	FILE *erpfile = fopen("Outputs/VentERP_TPORd_DOMINATE_ENDO.dat","w+");
	#endif


	//TP06* cell = new TP06(MCELL);
	//cell->init(MCELL);



	// apply user configuration about dt
	cell->setDt(dt);

	double time = 0;
	int step = 0;
	double oldV, oldDvdt, repo20, repo25, repo50, repo80, repo90;
	bool peakfound;
	for(time = 0.0, step = 0; time <= numS1*BCL; time += dt, step++)
	{
		// 1. Apply stimulus according to user configuration
		if(time - floor(time/BCL)*BCL >= stimStart && 
		   time - floor(time/BCL)*BCL < stimStart + stimDuration)	
		{
			cell->setIstim(stimStrength);
		}
		else
		{
			cell->setIstim(0.0);
		}

		// 2. update all states and currents
		cell->update();


		// 4. Statistics of single cell, including APD20, APD50, APD90, dVdt_max, resting potential, overshoot, amplitude		
		if (floor((time-stimStart)/BCL) == numS1 - 2 && floor((time+dt-stimStart)/BCL) == numS1 - 1)
		{
			rest_potential = cell->getV();
			oldDvdt = cell->getAbsDvdt();// record this to help detect the maximum dv/dt (UpstrokeVelocity_max)
			oldV = cell->getV();// record this to help detect whether the peak potential (overshoot)
			peakfound = false; // record whether peak potential (overshoot) found or not

			// 20190420 test
			dvdt_max = -1;
		}
		if (floor((time-stimStart)/BCL) == numS1 - 1)
		{
			if(cell->getV() > oldV)
			{
				// oldV = cell->getV();
				// dvdt_max = (cell->getAbsDvdt() > oldDvdt? cell->getAbsDvdt() : dvdt_max);
				// oldDvdt = cell->getAbsDvdt();

				// 20190420 test
				oldV = cell->getV();
				// dvdt_max = (cell->getAbsDvdt() > oldDvdt? cell->getAbsDvdt() : dvdt_max);
				// oldDvdt = cell->getAbsDvdt();


				// dvdt_max = (cell->getAbsDvdt() > dvdt_max? cell->getAbsDvdt() : dvdt_max);



				if( cell->getAbsDvdt() > dvdt_max )
				{
					dvdt_max = cell->getAbsDvdt();
					// cout << cell->getINa() << " " << cell->getIstim() << " " << cell->getICaL() << " " << cell->getINaK() << endl; 
				}


			}	
			else if(cell->getV() <= oldV && !peakfound) // peak not found yet
			{
				peakfound = true;
				overshoot = oldV;
				amplitude = overshoot - rest_potential; // should always be a positive value
				repo20 = overshoot - 0.20*amplitude;
				repo25 = overshoot - 0.25*amplitude;
				repo50 = overshoot - 0.50*amplitude;
				repo80 = overshoot - 0.80*amplitude;
				repo90 = overshoot - 0.90*amplitude;
				oldV = cell->getV();

			}
			else if(cell->getV() <= oldV && peakfound) // peak already found
			{
				if(oldV >= repo20 && cell->getV() <= repo20)
					apd20 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo25 && cell->getV() <= repo25)
					apd25 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo50 && cell->getV() <= repo50)
					apd50 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo80 && cell->getV() <= repo80)
					apd80 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo90 && cell->getV() <= repo90)
					apd90 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				oldV = cell->getV();
			}
		}

		// 5. output apd statistics to screen
		if (time + dt >= numS1*BCL)
		{	// APD20, APD50, APD90, ,  , , 
			printf("Resting membrane potential = %.5f mV\n",rest_potential); // unit: mV
			printf("Maximum dV/dt (UpstrokeVelocity_max) = %.5f mV/ms\n",dvdt_max); // unit: mV/ms
			printf("Overshoot = %.5f mV\n",overshoot); // unit: mV	
			printf("Amplitude = %.5f mV\n",amplitude); // unit: mV
			printf("AP20 = %.5f ms\n",repo20); // unit: ms
			printf("APD20 = %.5f ms\n",apd20); // unit: ms		
			printf("APD25 = %.5f ms\n",apd25); // unit: ms		
			printf("APD50 = %.5f ms\n",apd50); // unit: ms		
			printf("APD80 = %.5f ms\n",apd80); // unit: ms	
			printf("AP90 = %.5f ms\n",repo90); // unit: ms	
			printf("APD90 = %.5f ms\n",apd90); // unit: ms
		}
	}

	// S2 start.
	double stimS2start;
	bool erpfound = false;
	for(double DI = 0; DI >= -EPSILON; DI = DI + 1)
	{
		// re-initialize for every iteration.  
		// MUST BE RE-INITIALIZE AGAIN HERE!!!
		#ifdef EPI1
		//initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_EPI_antibody60.dat","r");
		initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_EPI.dat","r");
		cell->readinAllStates(initfile);
		#endif

		#ifdef ENDO1
		//initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_ENDO_antibody60.dat","r");
		initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_ENDO.dat","r");
		cell->readinAllStates(initfile);
		#endif

		#ifdef MCELL1
		//initfile = fopen("SingleCell/TP06InitialValues_DOMINATE_MCELL_antibody60.dat","r");
		initfile = fopen("SingleCell/TPORdInitialValues_DOMINATE_MCELL.dat","r");
		cell->readinAllStates(initfile);
		#endif

		stimS2start = stimStart + apd90 + DI;//start from APD90 to search for the refractory point.

		// write file                            ******可以选择不输出s1s2文件,要输出s1说文件的话需要手动更改文件名
		FILE *s1s2file;
		char s1s2filename[200];
		sprintf(s1s2filename, "Outputs/TPORd_ENDO_S1S2@%.1f.dat", fabs(stimS2start));
		s1s2file = fopen(s1s2filename,"w");
		

		// start a S2 at certain DI	
		double s2maxV = -1000;	
		for(time = 0, step = 0; time <= BCL; time += dt,step++)
		{

			erpfound = false;
			cell->setIstim(0);
			// stimulation applying or not
		    if(time >= stimStart && time < stimStart + stimDuration) // a normal but last S1
			{
				cell->setIstim(stimStrength);
			}
			else if(time >= stimS2start && 
				time < stimS2start + stimDuration) // S2 stim
			{
				cell->setIstim(stimStrength);
			}

			cell->update();
			//**************************************************控制是否需要输出文件，我们不需要膜电压的文件，我们只要ERP的值
			fprintf(s1s2file,"%4.10f\t",time);
			fprintf(s1s2file,"%4.10f\t",cell->getV());
			fprintf(s1s2file,"\n");

			if(time >= stimS2start)
			{
				if(cell->getV() > s2maxV) 
					s2maxV = cell->getV();// update max volt in s2
				if (cell->getV() >= repo20 && !erpfound)
				{
					erpfound = true;
					break;
				}
			}

		}

		fclose(s1s2file);



		if(erpfound)
		{
			cout << "ERP found! ERP = " << stimS2start - stimStart << " " << s2maxV<< endl ;
			fprintf(erpfile, "ERP found! ERP = %.10f   s2maxV= %.10f ",stimS2start - stimStart,s2maxV);
			break;
		}
		else
		{
			printf("S1S2 %.1f (ms) finished successfully. s2maxV = %.1f mV.\n", stimS2start - stimStart, s2maxV); 
		}

	} // S2 finish at certain BCL
	if(!erpfound) printf("No ERP found.\n");
	printf("All done!\n");


	fclose(initfile);
	//fclose(datafile);
	return 0;
}







