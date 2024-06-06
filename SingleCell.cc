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

int main(int argc, char *argv[])
{
	// --------user configuration list--------
	// All you need to do is put your single cell model into SingleCell folder
	// and modify following user-dependent parameters.
	
	#if defined(ATRIA) || defined(VENTRICLE)
	double BCL = 1000;//167; // ms
	double numS1 = 50; // 50 for short states. 49 for long states.
	double dt = 0.02; //ms
	double stopTime = numS1*BCL; //ms
	double stimStrength = -52;//-12;//-0.6e-3; // pA/pF   -6.0pA/pF(-0.6nA) for RatAtrial // -12.5pA/pF for CaMKII; -8.78pA for neonatalRatAtrial
	double stimDuration = 1;//3;   // ms
	double stimStart = 50.0;   // 20 ms  // indicates the time point of beginning stimulus in a cycle

	//****************************************我们的SingleCell中有两种单细胞模型，分别是TP06和TPORd，这里用的是TPORd。**********************
	//typedef TP06 CellType;
	//typedef TP06* CellPointer;
	typedef TPORd CellType;
	typedef TPORd* CellPointer;

	// statistics for single cell
	double apd20;
	double apd25;
	double apd50;
	double apd75;
	double apd80;
	double apd90;
	double dvdt_max;
	double rest_potential;
	double amplitude;
	double overshoot;
	#endif


	// --------start simulation--------
	//**********************手动修改细胞类型		2.需要手动修改输出文件的文件名


	// note that constructor contains the initializer
	CellPointer cell = new CellType(MCELL); // LVEPI for rats! EPI for human!                **********!!!!!!!MCELL要改成输出倒数第二个周期！！***
	// cell->setCORM2(30);
	// CellPointer cell = new CellType(EPI);
	// CellPointer cell = new CellType();
	#ifdef VENTRICLE
	FILE *datafile = fopen("Outputs/VentriSingleCellResults_TPORd_HOMO_MCELL.dat","w+");
	FILE *apdfile = fopen("Outputs/VentriSingleCellStats_TPORd_HOMO_MCELL.dat","w+");
	//FILE *initfile = fopen("SingleCell/TP06InitialValues_EPI.dat","w+");
	#endif

	
	// uncomment following two lines to read in initial values (this is because the original init values is not stable yet)
	// if the initfile is not available, run the initialization.cc first
	// FILE *initfile = fopen("SingleCell/TP06InitialValues_SQT6_EPI.dat","r");
	// cell->readinAllStates(initfile);


	// apply user configuration about dt
	cell->setDt(dt);

	double time = 0;
	double t_maxdvdt;
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
		if(step%(int(0.02/dt)) == 0 && floor(time/BCL) >= numS1 - 1) // if only focus on the last cycle              如果只输出最后一个周期
		//if(step%(int(0.1/dt)) == 0) // 5*dt = 1ms once
		//if(step%(int(0.02/dt)) == 0 && floor(time/BCL) >= numS1 - 2 && floor(time/BCL) < numS1 - 1)                //输出倒数第二个周期
		{
			 //fprintf(datafile,"%4.10f\t", time);                                                                   //输出所有周期
			fprintf(datafile,"%4.10f\t", fmod(time,BCL)); // if only focus on the last cycle*************************如果只输出最后一个周期
			// fprintf(datafile,"%4.10f\t", time - (numS1 - 2)*BCL ); // if focus on the last two

			// for alternans
			// fprintf(datafile,"%4.10f\t", cell->getV()); // unit: mV
			// fprintf(datafile,"%4.10f\t", cell->getNai()); // unit: mV
			// fprintf(datafile,"%4.10f\t", cell->getKi()); // unit: mV
			// fprintf(datafile,"%4.10f\t", cell->getCai()); // unit: mV

			// fprintf(datafile,"%4.10f\t", cell->getCaSR()); // unit: 
			// fprintf(datafile,"%4.10f\t", cell->getINaK()); // unit: 
			// fprintf(datafile,"%4.10f\t", cell->getISERCA()); // unit: 
			// fprintf(datafile,"%4.10f\t", cell->getIRyR()); // unit: 
			// fprintf(datafile,"%4.10f\t", cell->getISR()); // unit: 

			// for the affected current
			fprintf(datafile,"%4.10f\t", cell->getV()); // unit: mV
			// fprintf(datafile,"%4.10f\t", cell->getINa()); // unit: 
			fprintf(datafile,"%4.10f\t", cell->getICaL()); // unit: 
			// fprintf(datafile,"%4.10f\t", cell->getIK1()); // unit: 
			fprintf(datafile,"\n");
		}

		// 4. Statistics of single cell, including APD20, APD50, APD90, dvdt_max, resting potential, overshoot, amplitude
		double oldV, oldDvdt, repo20, repo25, repo50, repo75, repo80, repo90;
		bool peakfound;
		// time-stimStart is used here.
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
				repo75 = overshoot - 0.75*amplitude;
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
				else if(oldV >= repo75 && cell->getV() <= repo75)
					apd75 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo80 && cell->getV() <= repo80)
					apd80 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				else if(oldV >= repo90 && cell->getV() <= repo90)
					apd90 = time - floor(time/BCL)*BCL - stimStart - stimDuration; // note that apd calculate from the stimStart.
				oldV = cell->getV();
			}
		}

		// 5. output apd statistics to file
		if (time + dt >= stopTime)
		{	// APD20, APD50, APD90, ,  , , 
			fprintf(apdfile,"Resting membrane potential = %.5f mV\n",rest_potential); // unit: mV
			fprintf(apdfile,"Maximum dV/dt (UpstrokeVelocity_max) = %.5f mV/ms\n",dvdt_max); // unit: mV/ms
			fprintf(apdfile,"Overshoot = %.5f mV\n",overshoot); // unit: mV	
			fprintf(apdfile,"Amplitude = %.5f mV\n",amplitude); // unit: mV
			fprintf(apdfile,"APD20 = %.5f ms\n",apd20); // unit: ms					
			fprintf(apdfile,"APD25 = %.5f ms\n",apd25); // unit: ms	
			fprintf(apdfile,"APD50 = %.5f ms\n",apd50); // unit: ms	
			fprintf(apdfile,"APD75 = %.5f ms\n",apd75); // unit: ms		
			fprintf(apdfile,"APD80 = %.5f ms\n",apd80); // unit: ms		
			fprintf(apdfile,"APD90 = %.5f ms\n",apd90); // unit: ms
		}

		// 6. output apd statistics to screen
		if (time + dt >= stopTime)
		{	// APD20, APD50, APD90, ,  , , 
			printf("Resting membrane potential = %.5f mV\n",rest_potential); // unit: mV
			printf("Maximum dV/dt (UpstrokeVelocity_max) = %.5f mV/ms\n",dvdt_max); // unit: mV/ms
			printf("Overshoot = %.5f mV\n",overshoot); // unit: mV	
			printf("Amplitude = %.5f mV\n",amplitude); // unit: mV
			printf("APD20 = %.5f ms\n",apd20); // unit: ms	
			printf("AP20 = %.5f ms\n",repo20); // unit: ms	
			printf("APD25 = %.5f ms\n",apd25); // unit: ms	
			printf("AP25 = %.5f ms\n",repo25); // unit: ms		
			printf("APD50 = %.5f ms\n",apd50); // unit: ms	
			printf("APD75 = %.5f ms\n",apd75); // unit: ms		
			printf("APD80 = %.5f ms\n",apd80); // unit: ms		
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

			// cell->outputAllStates(co_initfile); // to get the alternans initial states.
		}
	}
	fclose(datafile);
	fclose(apdfile);
	// fclose(co_initfile);
	return 0;
}




