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
 * Last update : 06-10-2023
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


//#define wildtype
//#define hete
#define homo

//#define epi1
//#define endo1
#define mcell1

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	// rate depedency
	bool finish = false;
	double BCL = 4000;
	double numS1 = 50;//49,50

#if defined (wildtype) && defined(epi1)
	FILE* rcfile = fopen("Outputs/RestCurve_wild-type_epi_S1S2.dat", "w+");
#endif // define wild-type && epi1
#if defined (wildtype) && defined(endo1)
	FILE* rcfile = fopen("Outputs/RestCurve_wild-type_endo_S1S2.dat", "w+");
#endif // define wild-type && endo1
#if defined (wildtype) &&defined (mcell1)
	FILE* rcfile = fopen("Outputs/RestCurve_wild-type_mcell_S1S2.dat", "w+");
#endif // define wild-type && endo1

#if defined (hete) &&defined( epi1)
	FILE* rcfile = fopen("Outputs/RestCurve_hete_epi_S1S2.dat", "w+");
#endif // define hete && epi1
#if defined (hete) &&defined( endo1)
	FILE* rcfile = fopen("Outputs/RestCurve_hete_endo_S1S2.dat", "w+");
#endif // define hete && endo1
#if defined (hete) &&defined(mcell1)
	FILE* rcfile = fopen("Outputs/RestCurve_hete_mcell_S1S2.dat", "w+");
#endif // define hete &&

#if defined(homo)&&defined(epi1)
	FILE* rcfile = fopen("Outputs/RestCurve_homo_epi_S1S2.dat", "w+");
#endif // defined(homo)&&defined(epi1)
#if defined (homo) &&defined( endo1)
	FILE* rcfile = fopen("Outputs/RestCurve_homo_endo_S1S2.dat", "w+");
#endif // define hete && endo1
#if defined (homo) &&defined(mcell1)
	FILE* rcfile = fopen("Outputs/RestCurve_homo_mcell_S1S2.dat", "w+");
#endif 

	while(!finish)
	{
		#if defined(ATRIA) || defined(VENTRICLE)
		double dt = 0.02; //ms
		double stopTime = numS1*BCL; //ms
		double stimStrength = -52.0; // pA/pF   -6.0pA/pF(-0.6nA) for RatAtrial // -12.5pA/pF for CaMKII; -8.78pA for neonatalRatAtrial
		double stimDuration = 1;   // ms
		double stimStart = 0.0;   // ms  // do not set this value. it will lead to incorrect calculation
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


		
		// uncomment following two lines to read in initial values (this is because the original init values is not stable yet)
		// if the initfile is not available, run the initialization.cc first
		// FILE *initfile = fopen("SingleCell/TP06InitialValues_CON_EPI.dat","r");
		// cell->readinAllStates(initfile);
#if defined (wildtype) && defined(epi1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_WT_EPI.dat", "r");
		
#endif // define wild-type && epi1
#if defined (wildtype) && defined(endo1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_WT_ENDO.dat", "r");
	
#endif // define wild-type && endo1
#if defined (wildtype) &&defined (mcell1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_WT_MCELL.dat", "r");
		
#endif // define wild-type && endo1

#if defined (hete) &&defined( epi1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HETE_EPI.dat", "r");
		
#endif // define hete && epi1
#if defined (hete) &&defined( endo1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HETE_ENDO.dat", "r");
		
#endif // define hete && endo1
#if defined (hete) &&defined(mcell1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HETE_MCELL.dat", "r");
	
#endif // define hete &&

#if defined(homo)&&defined(epi1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HOMO_EPI.dat", "r");
	
#endif // defined(homo)&&defined(epi1)
#if defined (homo) &&defined( endo1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HOMO_ENDO.dat", "r");
		
#endif // define hete && endo1
#if defined (homo) &&defined(mcell1)
		FILE* initfile = fopen("SingleCell/TP06InitialValues_HOMO_MCELL.dat", "r");
	
#endif 

		cell->readinAllStates(initfile);
		// apply user configuration about dt
		cell->setDt(dt);

		double time = 0;
		int step = 0;
		for(time = 0.0, step = 0; time < stopTime; time += dt, step++)
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
			double oldV, oldDvdt, repo20, repo25, repo50, repo80, repo90;
			bool peakfound;
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
					oldV = cell->getV();

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

			// 5. output apd statistics to file
			if (time + dt >= stopTime)
			{	// APD20, APD50, APD90
				fprintf(rcfile,"%.5f\t",BCL); // unit: ms		
				fprintf(rcfile,"%.5f\n",apd90); // unit: ms
			}

			// 6. output apd statistics to screen
			if (time + dt >= stopTime)
			{	// APD20, APD50, APD90, ,  , , 
				// printf("Resting membrane potential = %.5f mV\n",rest_potential); // unit: mV
				// printf("Maximum dV/dt (UpstrokeVelocity_max) = %.5f mV/ms\n",dvdt_max); // unit: mV/ms
				// printf("Overshoot = %.5f mV\n",overshoot); // unit: mV	
				// printf("Amplitude = %.5f mV\n",amplitude); // unit: mV
				// printf("APD20 = %.5f ms\n",apd20); // unit: ms		
				// printf("APD25 = %.5f ms\n",apd25); // unit: ms		
				// printf("APD50 = %.5f ms\n",apd50); // unit: ms		
				// printf("APD80 = %.5f ms\n",apd80); // unit: ms	
				printf("BCL = %.5f ms\t",BCL);	
				printf("APD90 = %.5f ms\n",apd90); // unit: ms

				/*
				cout << cell->getV()        << endl;//       {return Volt;}
				cout << cell->getM()        << endl;//        {return m;}
				cout << cell->getH()        << endl;//        {return h;}
				cout << cell->getJ()        << endl;//        {return j;}
				cout << cell->getD()        << endl;//        {return d;}
				cout << cell->getF11()      << endl;//        {return f11;}
				cout << cell->getF12()      << endl;//        {return f12;}
				cout << cell->getFCa()      << endl;//        {return Cainact;}
				cout << cell->getR()        << endl;//            {return r;}
				cout << cell->getS()        << endl;//            {return s;}
				cout << cell->getSslow()    << endl;//                {return sslow;}
				cout << cell->getRss()      << endl;//              {return rss;}
				cout << cell->getSss()      << endl;//              {return sss;}
				cout << cell->getY()        << endl;//            {return y;}
				cout << cell->getNai()      << endl;//              {return Nai;}
				cout << cell->getKi()       << endl;//             {return Ki;}
				cout << cell->getCai()      << endl;//              {return Cai;}
				cout << cell->getCaNSR()    << endl;//                {return CaNSR;}
				cout << cell->getCaSS()     << endl;//               {return CaSS;}
				cout << cell->getCaJSR()    << endl;//                {return CaJSR;}
				cout << cell->getPC1()      << endl;//              {return PC1;}
				cout << cell->getPo1()      << endl;//              {return Po1;}
				cout << cell->getPo2()      << endl;//              {return Po2;}
				cout << cell->getPC2()      << endl;//              {return PC2;}
				cout << cell->getLTRPNCa()  << endl;//                  {return LTRPNCa;}
				cout << cell->getHTRPNCa()  << endl;//                  {return HTRPNCa;}
				*/
			}
		}
		




		// the dynamic protocol: prepare the next round.
		if(BCL > 4000)
	    {
			BCL = BCL - 1000;
			// numS1 = 50;
	    }
	  	else if(BCL > 3000)
	    {
	      	BCL = BCL - 1000;
	      	// numS1 = 50;
	    }
	  	else if(BCL > 2000)
	    {
	      	BCL = BCL - 200;
	      	// numS1 = 50;
	    }
	  	else if(BCL > 1000)
	    {
	      	BCL = BCL - 100;
	      	// numS1 = 50;
	    }
	  	else if(BCL > 500)
	    {	     
	      	BCL = BCL - 50;
	      	// numS1 = 50;
	    }
	  	else if(BCL > 400)
	    {
	      	BCL = BCL - 10;
	      	// numS1 = 50;
	    }
	  	else if(BCL > 300)
	    {
	      	BCL = BCL - 5;
	      	// numS1 = 50;
	    }
	  	else if(BCL > 250)
	    {
	      	BCL = BCL - 5;
	      	// numS1 = 50;
	    }
	  	else
	    {
	    	finish = true;
	      	printf("Restitution protocol finished.\n");
	    }
	}
	fclose(rcfile);
	return 0;
}




