/*
 * General Code Structure (GCS) for single cell simulation
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


// #include "SingleCell/NeonatalRatAtria.cc"
// #include "SingleCell/RatSAN.cc"
// #include "SingleCell/GrandiCaMK.cc"
// #include "SingleCell/HumanVentricle.cc"
// #include "SingleCell/RatLV.h"
#include "SingleCell/TP06.h"

using namespace std;

// ONLY ONE LINE CAN BE KEPT HERE
// #define SAN
// #define ATRIA
#define VENTRICLE


#define wildtype
//#define hete
//#define homo

#define epi1
//#define endo1
//#define mcell1



int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	// rate depedency
	double BCL = 1000;
	double numS1 = 2;
/*
	#ifdef wildtype
	FILE *rcfile = fopen("Outputs/RestCurve_wild-type_S1S2.dat","w+");
	#endif
	#ifdef hete
	FILE *rcfile = fopen("Outputs/RestCurve_hete_S1S2.dat","w+");
	#endif
	#ifdef homo
	FILE *rcfile = fopen("Outputs/RestCurve_homo_S1S2.dat","w+");
	#endif
*/

	#if defined(ATRIA) || defined(VENTRICLE)
	double dt = 0.02; //ms
	double s1stopTime = numS1*BCL; //ms
	double stimStrength = -52.0; // pA/pF   -6.0pA/pF(-0.6nA) for RatAtrial // -12.5pA/pF for CaMKII; -8.78pA for neonatalRatAtrial
	double stimDuration = 1;   // ms
	double stimStart = 50.0;   // ms  // indicates the time point of beginning stimulus in a cycle
	// typedef GrandiCaMKII CellType;
	// typedef GrandiCaMKII* CellPointer;
	// typedef ORdHumanVentricle CellType;
	// typedef ORdHumanVentricle* CellPointer;
	// typedef NeonatalRatAtria CellType;
	// typedef NeonatalRatAtria* CellPointer;
	typedef TP06 CellType;
	typedef TP06* CellPointer;

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

	double oldV, oldDvdt, repo20, repo25, repo50, repo80, repo90;
	bool peakfound;
	#endif



	// --------start simulation--------



	// note that constructor contains the initializer
	#ifdef epi1
		CellPointer cell = new CellType(EPI);
	#endif // epi1
	#ifdef endo1
		CellPointer cell = new CellType(ENDO);
	#endif // endo1
	#ifdef mcell1
		CellPointer cell = new CellType(MCELL);
	#endif // mcell1



	//CellPointer cell = new CellType(EPI);
	// CellPointer cell = new CellType();


	// uncomment following two lines to read in initial values (this is because the original init values is not stable yet)
	// if the initfile is not available, run the initialization.cc first
#if defined (wildtype) && defined(epi1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_WT_EPI.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_wild-type_epi_S1S2.dat", "w+");
#endif // define wild-type && epi1
#if defined (wildtype) && defined(endo1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_WT_ENDO.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_wild-type_endo_S1S2.dat", "w+");
#endif // define wild-type && endo1
#if defined (wildtype) &&defined (mcell1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_WT_MCELL.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_wild-type_mcell_S1S2.dat", "w+");
#endif // define wild-type && endo1

#if defined (hete) &&defined( epi1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HETE_EPI.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_hete_epi_S1S2.dat", "w+");
#endif // define hete && epi1
#if defined (hete) &&defined( endo1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HETE_ENDO.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_hete_endo_S1S2.dat", "w+");
#endif // define hete && endo1
#if defined (hete) &&defined(mcell1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HETE_MCELL.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_hete_mcell_S1S2.dat", "w+");
#endif // define hete &&

#if defined(homo)&&defined(epi1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HOMO_EPI.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_homo_epi_S1S2.dat", "w+");
#endif // defined(homo)&&defined(epi1)
#if defined (homo) &&defined( endo1)
			FILE* initfile = fopen("SingleCell/TP06InitialValues_HOMO_ENDO.dat", "r");
			FILE* rcfile = fopen("Outputs/RestCurve_homo_endo_S1S2.dat", "w+");
#endif // define hete && endo1
#if defined (homo) &&defined(mcell1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HOMO_MCELL.dat", "r");
		FILE* rcfile = fopen("Outputs/RestCurve_homo_mcell_S1S2.dat", "w+");
#endif // define hete &&


	//FILE *initfile = fopen("SingleCell/TP06InitialValues_CON_EPI.dat","r");
	cell->readinAllStates(initfile);


	// apply user configuration about dt
	cell->setDt(dt);

	double time = 0;
	int step = 0;
	for(time = 0.0, step = 0; time < s1stopTime; time += dt, step++)
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

		// 3. Output file, write file each 1ms 


		// 4. Statistics of single cell, including APD20, APD50, APD90, dVdt_max, resting potential, overshoot, amplitude
		if (floor((time-stimStart)/BCL) == numS1 - 2 && floor(((time-stimStart)+dt)/BCL) == numS1 - 1)
		{
			rest_potential = cell->getV();
			oldDvdt = cell->getAbsDvdt();// record this to help detect the maximum dv/dt (UpstrokeVelocity_max)
			oldV = cell->getV();// record this to help detect whether the peak potential (overshoot)
			peakfound = false; // record whether peak potential (overshoot) found or not
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



				if(cell->getAbsDvdt() > dvdt_max )
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
				// cout << "s1 overshoot = " << overshoot << endl;

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

	}// s1 finished




	// output s1 statistics to screen
	printf("S1 stimuli finished! APD90 = %.5f ms\n",apd90); // unit: ms

	// save s1 results for s2 initialization
	FILE *s2_initfile = fopen("Outputs/S2_INITFILE.dat","w+");
	cell->outputAllStates(s2_initfile);
	// fclose(s2_initfile);

	FILE *datafile = fopen("Outputs/s2data.dat","w+");

	// loop s2
	double stimS2start = 0;
	double DI = 2000;
	while(DI > 50)
	{
		// init
		rewind(s2_initfile);
		cell->readinAllStates(s2_initfile);

		
		// the dynamic protocol: prepare the next round.
		if(DI > 4000) 
		{
			DI = DI - 1000;
		}
		else if(DI > 3000)
		{
		  	DI = DI - 1000;
		}
		else if(DI > 2000)
		{
		  	DI = DI - 1000;
		}
		else if(DI > 1000)
		{
		  	DI = DI - 1000;
		}
		else if(DI > 500)
		{	     
		  	DI = DI - 250;
		}
		else if(DI > 400)
		{
		  	DI = DI - 50;
		}
		else if(DI > 300)
		{
		  	DI = DI - 10;
		}
		else if(DI > 250)
		{
		  	DI = DI - 5;
		}
		else if(DI > 50)
		{
		  	DI = DI - 1;
		}

		stimS2start = stimStart + stimDuration + apd90 + DI;
		
		double s2repo90, s2apd90;
		rest_potential = cell->getV();
		oldV = cell->getV();
		overshoot = repo90;
		// cout << "1. overshoot " << overshoot << endl;
		peakfound = false;

		

		for(time = 0.0, step = 0; time < apd90 + DI + BCL; time += dt, step++)
		{
			if((time >= stimStart && time < stimStart + stimDuration) ||     // s1
			   (time >= stimS2start && time < stimS2start + stimDuration))  // s2
			{
				cell->setIstim(stimStrength);
			}
			else
			{
				cell->setIstim(0.0);
			}

			// 2. update all states and currents
			cell->update();

			if(DI >= 749 && DI < 751 && step%(int(1/dt)) == 0) // 5*dt = 1ms once
			{
				fprintf(datafile,"%4.10f\t", time);
				fprintf(datafile,"%4.10f\n", cell->getV()); // unit: mV
			}
			

			// now is s2 period
			if(time >= stimS2start)
			{
				if(cell->getV() > oldV)
				{
					oldV = cell->getV();
					// cout << "1. oldV = " << oldV << endl;
				}

				else if(cell->getV() <= oldV && !peakfound) // peak not found yet
				{
					peakfound = true;					
					overshoot = oldV;
					// cout << "2. oldV = " << oldV << endl;
					// cout << "overshoot = " << overshoot << endl;
					amplitude = overshoot - rest_potential; // should always be a positive value
					s2repo90 = overshoot - 0.90*amplitude;
					oldV = cell->getV();


				}

				else if(cell->getV() <= oldV && peakfound) // peak already found
				{
					if(oldV >= s2repo90 && cell->getV() <= s2repo90)
					{
						s2apd90 = time - stimS2start - stimDuration; // note that apd calculate from the end of the S2 stimDuration.
						cout << "DI = " << DI << "	APD90 = " << s2apd90 << endl;
						fprintf(rcfile,"%.5f\t",DI); // unit: ms
						fprintf(rcfile,"%.5f\n",s2apd90); // unit: ms
						break; // for next DI.
					}
					oldV = cell->getV();
				}
			} // s2 finished.

		} // certain DI finished.


	} // DI loop finished.


	fclose(s2_initfile);
	fclose(rcfile);
	fclose(datafile);
	return 0;
}




