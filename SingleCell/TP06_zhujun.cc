#ifndef _TP06_CC_
#define _TP06_CC_
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <typeinfo>
#include "omp.h"
#include <ctime>
#include "Cell.cc"
#include "TP06.h"

using namespace std;


// NOTE: ONLY ONE OF THE FOLLOWING THREE LINES CAN BE KEPT HERE!
// #define BASELINE // original TP06 model
#define CON      // Control for casestudy
// #define CORM


// Description can be found in TP06.h

TP06::TP06(CellType ct)
{
	init(ct);
}

void TP06::init(CellType ct)
{
	//Initial values of state variables
	svolt = -86.2;
	Cai = 0.00007;
	CaSR = 1.3;
	CaSS = 0.00007;
	Nai = 7.67;
	Ki = 138.3;
	sm = 0.;
	sh = 0.75;
	sj = 0.75; 
	mL = 0;
	hL = 1;
	hLp = 1;
	sxr1 = 0.;
	sxr2 = 1.;
	sxs = 0.;
	sr = 0.;
	ss = 1.;
	sd = 0.;
	sf = 1.;
	sf2 =1.;
	sfcass = 1.;                                        
	sRR = 1.;  
	sOO = 0.;

	ctype = ct;
}

void TP06::outputAllStates(FILE *datafile)
{
	fprintf(datafile,"%4.10f\n", svolt   );
    fprintf(datafile,"%4.10f\n", Cai     );
    fprintf(datafile,"%4.10f\n", CaSR    );
    fprintf(datafile,"%4.10f\n", CaSS    );
    fprintf(datafile,"%4.10f\n", Nai     );
    fprintf(datafile,"%4.10f\n", Ki      );
    fprintf(datafile,"%4.10f\n", sm      );
    fprintf(datafile,"%4.10f\n", sh      );
    fprintf(datafile,"%4.10f\n", sj      );
    fprintf(datafile,"%4.10f\n", mL      );
    fprintf(datafile,"%4.10f\n", hL      );
    fprintf(datafile,"%4.10f\n", hLp     );   
    fprintf(datafile,"%4.10f\n", sxr1    );
    fprintf(datafile,"%4.10f\n", sxr2    );
    fprintf(datafile,"%4.10f\n", sxs     );
    fprintf(datafile,"%4.10f\n", sr      );
    fprintf(datafile,"%4.10f\n", ss      );
    fprintf(datafile,"%4.10f\n", sd      );
    fprintf(datafile,"%4.10f\n", sf      );
    fprintf(datafile,"%4.10f\n", sf2     );
    fprintf(datafile,"%4.10f\n", sfcass  );
    fprintf(datafile,"%4.10f\n", sRR     );
    fprintf(datafile,"%4.10f\n", sOO     );
}

void TP06::readinAllStates(FILE *datafile)
{
	double value;
	fscanf(datafile,"%lf", &value); 
	svolt  = value;
	fscanf(datafile,"%lf", &value); 
	Cai    = value;
	fscanf(datafile,"%lf", &value); 
	CaSR   = value;
	fscanf(datafile,"%lf", &value); 
	CaSS   = value;
	fscanf(datafile,"%lf", &value); 
	Nai    = value;
	fscanf(datafile,"%lf", &value); 
	Ki     = value;
	fscanf(datafile,"%lf", &value); 
	sm     = value;
	fscanf(datafile,"%lf", &value); 
	sh     = value;
	fscanf(datafile,"%lf", &value); 
	sj     = value;
	fscanf(datafile,"%lf", &value); 
	mL     = value;
	fscanf(datafile,"%lf", &value); 
	hL     = value;
	fscanf(datafile,"%lf", &value); 
	hLp    = value;
	fscanf(datafile,"%lf", &value); 
	sxr1   = value;
	fscanf(datafile,"%lf", &value); 
	sxr2   = value;
	fscanf(datafile,"%lf", &value); 
	sxs    = value;
	fscanf(datafile,"%lf", &value); 
	sr     = value;
	fscanf(datafile,"%lf", &value); 
	ss     = value;
	fscanf(datafile,"%lf", &value); 
	sd     = value;
	fscanf(datafile,"%lf", &value); 
	sf     = value;
	fscanf(datafile,"%lf", &value); 
	sf2    = value;
	fscanf(datafile,"%lf", &value); 
	sfcass = value;
	fscanf(datafile,"%lf", &value); 
	sRR    = value;
	fscanf(datafile,"%lf", &value); 
	sOO    = value;
}

void TP06::update()
{
	//void Step(Variables *V,double dt,char *despath,double *tt,int step,double Istim)
	//#define v(array_pointer,i,j) (*(V->array_pointer+i*V->NJ +j))

	// those following variables can be simply regarded as temp values, and no need to be put in header file
	/*-----------------------------------------------------------------------------
	ELECTROPHYSIOLOGICAL PARAMETERS:可以把这些都放在update里
	-----------------------------------------------------------------------------*/
	double const ATP=5;
	double const ADP=75;
	double const CO=15;//CO concentration μM
	double tK1= log10(CO*1e-6);
	double tKr= log10(CO*1e-6);
	double tCaL= log10(CO);
	double tNa= log10(CO*1e-6);

	double RK1= 0.9/(1.0+pow(-5.256/tK1,5.42))+0.1;
	double RKr= 0.73/(1.0+pow(-5.504/tKr,20.301))+0.265;
	double RCaL= 0.511/(1.0+pow(1.201/tCaL,-6.053))+0.4;
	double RNa= 1.2/(1.0+pow(-6.007/tNa,10.745))+0.0;

	//External concentrations
	double Ko = 5.4;
	// double Ko = 8.0;// 缺血
	double Cao = 2.0;
	double Nao = 140.0;

	//Intracellular volumes
	double Vc = 0.016404;
	double Vsr = 0.001094;
	double Vss = 0.00005468;

	//Calcium buffering dynamics
	double Bufc = 0.2;
	double Kbufc = 0.001;
	double Bufsr = 10.;
	double Kbufsr = 0.3;
	double Bufss = 0.4;
	double Kbufss = 0.00025;

	//Intracellular calcium flux dynamics
	double Vmaxup = 0.006375;
	double Kup = 0.00025;
	double Vrel = 0.102;//40.8;
	//double Vrel = 40.8;
	double k1_ = 0.15;
	double k2_ = 0.045;
	double k3 = 0.060;
	double k4 = 0.005;//0.000015;
	//double k4 = 0.000015;
	double EC = 1.5;
	double maxsr = 2.5;
	double minsr = 1.;
	double Vleak = 0.00036;
	double Vxfer = 0.0038;



	//Constants
	double R = 8314.472;
	double F = 96485.3415;
	double T = 310.0;
	double RTONF = (R*T)/F;

	//Cellular capacitance         
	double CAPACITANCE = 0.185;

	//Parameters for currents
	//Parameters for IKr
	double Gkr = 0.153;
	//double Gkr = 0.172;
		
	//Parameters for Iks
	double pKNa = 0.03;
	double Gks;
	if(ctype == MCELL)  
		Gks = 0.098;
	else // EPI & ENDO
		Gks = 0.392;
		//Gks = 0.441;
		
	//Parameters for Ik1
	double GK1 = 5.405;
	//Parameters for Ito
	double Gto;
	if(ctype == ENDO) 
		Gto = 0.073;
	else // EPI & MCELL
		Gto = 0.294;
	//Parameters for INa
	double GNa = 14.838;
	// double GNa = 14.838*0.8;//缺血	
	//Parameters for IbNa
	double GbNa = 0.00029;
	//Parameters for INaK
	double KmK = 1.0;
	double KmNa = 40.0;
	double knak = 2.724;
	//Parameters for ICaL
	#ifdef BASELINE
	double GCaL = 0.00003980;
	# endif
	#if defined (CON) || defined (CORM)
	// scaling factor 0.838 was applied to achieve same peak density as in original TP06 model
	double GCaL = 0.838*0.00003980; // Control for casestudy
	// double GCaL = 0.838*0.00003980*0.8; // 缺血
 	#endif 
	#ifdef SQT6
	double GCaL = 0.838*0.23*0.00003980; // SQT6 for casestudy
	#endif 
	//Parameters for IbCa
	double GbCa = 0.000592;
	//Parameters for INaCa
	double knaca = 1000;
	double KmNai = 87.5;
	double KmCa = 1.38;
	double ksat = 0.1;
	double n = 0.35;
	//Parameters for IpCa
	double GpCa = 0.1238;
	//double GpCa = 0.8666;
	
	double KpCa = 0.0005;
	//Parameters for IpK;
	double GpK = 0.0146;
	//double GpK = 0.00219;	




	// current	
	double k1;
	double k2;
	double kCaSR;

	// derivatives
	double dNai,dKi,dCai,dCaSR,dCaSS,dRR;

 	// reverse potential
	double Ek,Ena,Eks,Eca;

	double CaCSQN;
	double bjsr;
	double cjsr;
	double CaSSBuf;
	double bcss;
	double ccss;
	double CaBuf;
	double bc;
	double cc;
	double Ak1;
	double Bk1;
	double rec_iK1;
	double rec_ipK;
	double rec_iNaK;
	double AM;
	double BM;
	double AH_1;
	double BH_1;
	double AH_2;
	double BH_2;
	double AJ_1;
	double BJ_1;
	double AJ_2;
	double BJ_2;
	double M_INF;
	double H_INF;
	double J_INF;
	double TAU_M;
	double TAU_H;
	double TAU_J;
	double axr1;
	double bxr1;
	double axr2;
	double bxr2;
	double Xr1_INF;
	double Xr2_INF;
	double TAU_Xr1;
	double TAU_Xr2;
	double Axs;
	double Bxs;
	double Xs_INF;
	double TAU_Xs;
	double R_INF;
	double TAU_R;
	double S_INF;
	double TAU_S;
	double Ad;
	double Bd;
	double Cd;
	double Af;
	double Bf;
	double Cf;
	double Af2;
	double Bf2;
	double Cf2;
	double TAU_D;
	double D_INF;
	double TAU_F;
	double F_INF;
	double TAU_F2;
	double F2_INF;
	double TAU_FCaSS;
	double FCaSS_INF;

	
	static double inverseVcF2 = 1/(2*Vc*F);
	static double inverseVcF = 1./(Vc*F);
	static double inversevssF2 = 1/(2*Vss*F);

	 
	static char s[200];
	FILE *FF;

  
    
    //Needed to compute currents
    Ek = RTONF*(log((Ko/Ki)));
    Ena = RTONF*(log((Nao/Nai)));
    Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Eca = 0.5*RTONF*(log((Cao/Cai)));
    // Eca = 50; // case study voltage clamp
    Ak1 = 0.1/(1.+exp(0.06*(svolt-Ek-200)));
    Bk1 = (3.*exp(0.0002*(svolt-Ek+100))+
	 exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
    rec_iK1 = Ak1/(Ak1+Bk1);
    rec_iNaK = (1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
    rec_ipK = 1./(1.+exp((25-svolt)/5.98));


	double mLss=1.0/(1.0+exp((-(svolt+42.85))/5.264));
	double tm=1.0/(6.765*exp((svolt+11.64)/34.77)+8.552*exp(-(svolt+77.42)/5.955));
	double tmL=tm;
	mL=mLss-(mLss-mL)*exp(-dt/tmL);
	double hLss=1.0/(1.0+exp((svolt+87.61)/7.488));
	double thL=200.0;
	hL=hLss-(hLss-hL)*exp(-dt/thL);
	double hLssp=1.0/(1.0+exp((svolt+93.81)/7.488));
	double thLp=3.0*thL;
	hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
	double GNaL=0.0075;
	if (ctype == EPI)  //endo = 0, epi = 1, M = 2
	{
		GNaL*=0.6;
	}
	//double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
	double IKATP;
	double KMg, fM, KNa, fN, fT, Km, H, fATP;

	KMg = (0.65 / pow((Ko + 5), 0.5))*exp(-2 * 0.32/RTONF*svolt);
	fM = 1 / (1 + (3.1 / KMg));

	KNa = 25.9*exp(-0.35/RTONF*svolt);
	fN = 1 / (1 + pow((8 / KNa), 2));

	fT = pow(1.3, (T - 309.15) / 10);

	Km = 35.8 + 17.9*pow(ADP, 0.256);
	H = 1.3 + 0.74*exp(-0.09*ADP);
	fATP = 1 / (1 + pow((1000 * ATP / Km), H));

	IKATP = 3.495*pow(Ko / 5.4, 0.24)*0.91*fM*fN*fT*fATP*(svolt - Ek);

		
    //Compute currents
    #ifdef CON
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    //INaL=GNaL*(v-Ena)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
    INaL=GNaL*(svolt-Ena)*mL*hL;
    ICaL = GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
    (0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.);
    IKr =Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);    //debug  100% 50% 0%
    IK1 = GK1*rec_iK1*(svolt-Ek); // different from Lu
	#endif
 
    #ifdef CORM
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena)*0.5;
    //INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
    INaL=GNaL*(svolt-Ena)*mL*hL*2;
    ICaL = GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*(0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.)*0.6;
    IKr = Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek)*0.3;
    IK1 = GK1*rec_iK1*(svolt-Ek)*0.3;
	#endif


    Ito = Gto*sr*ss*(svolt-Ek);
    IKs = Gks*sxs*sxs*(svolt-Eks);     
    INaCa = knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
      (exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
       exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa = GpCa*Cai/(KpCa+Cai);
    IpK = GpK*rec_ipK*(svolt-Ek);
    IbNa = GbNa*(svolt-Ena);
    IbCa = GbCa*(svolt-Eca);


    //Determine total current   pA/pF
    dvdt = -(IKr +
      IKs   +
      IK1   +
      Ito   +
      INa   +
      INaL  +
      IbNa  +
      ICaL  +
      IbCa  +
      INaK  +
      INaCa +
      IpCa  +
      IpK   +
//       IKATP +
      Istim);

 
    //update concentrations    
    kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
    k1=k1_/kCaSR;
    k2=k2_*kCaSR;
    dRR=k4*(1-sRR)-k2*CaSS*sRR;
    sRR+=dt*dRR;
    sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);


    Irel=Vrel*sOO*(CaSR-CaSS);
    Ileak=Vleak*(CaSR-Cai);
    Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
    Ixfer=Vxfer*(CaSS-Cai);
    

    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=dt*(Iup-Irel-Ileak);
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
   

    CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
    dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
    bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
    ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
    CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;


    CaBuf=Bufc*Cai/(Cai+Kbufc);
    dCai=dt*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc=Bufc-CaBuf-dCai-Cai+Kbufc;
    cc=Kbufc*(CaBuf+dCai+Cai);
    Cai=(sqrt(bc*bc+4*cc)-bc)/2;
        
    
    dNai = -(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
    Nai += dt*dNai;
    
    dKi = -(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    Ki += dt*dKi;



    //compute steady state values and time constants 
    AM=1./(1.+exp((-60.-svolt)/5.));
    BM=0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
    TAU_M=AM*BM;
    M_INF=1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
    if (svolt>=-40.)
    {
		AH_1=0.; 
		BH_1=(0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
		TAU_H= 1.0/(AH_1+BH_1);
	}
    else
    {
		AH_2=(0.057*exp(-(svolt+80.)/6.8));
		BH_2=(2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
		TAU_H=1.0/(AH_2+BH_2);
    }
    H_INF=1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
    if(svolt>=-40.)
    {
		AJ_1=0.;      
		BJ_1=(0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
		TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
		AJ_2=(((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*
			exp(-0.04391*svolt))*(svolt+37.78)/
		      (1.+exp(0.311*(svolt+79.23))));    
		BJ_2=(0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
		TAU_J= 1.0/(AJ_2+BJ_2);
    }
    J_INF=H_INF;

    Xr1_INF=1./(1.+exp((-26.-svolt)/7.));
    axr1=450./(1.+exp((-45.-svolt)/10.));
    bxr1=6./(1.+exp((svolt-(-30.))/11.5));
    TAU_Xr1=axr1*bxr1;
    Xr2_INF=1./(1.+exp((svolt-(-88.))/24.));
    axr2=3./(1.+exp((-60.-svolt)/20.));
    bxr2=1.12/(1.+exp((svolt-60.)/20.));
    TAU_Xr2=axr2*bxr2;

    Xs_INF=1./(1.+exp((-5.-svolt)/14.));
    Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
    Bxs=(1./(1.+exp((svolt-35.)/15.)));
    TAU_Xs=Axs*Bxs+80;
    
	if(ctype == EPI)
	{
	    R_INF=1./(1.+exp((20-svolt)/6.));
	    S_INF=1./(1.+exp((svolt+20)/5.));
	    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
	    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
	}
	else if(ctype == ENDO) 
	{
	    R_INF=1./(1.+exp((20-svolt)/6.));
	    S_INF=1./(1.+exp((svolt+28)/5.));
	    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
	    TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
	}
	else // MCELL
	{
		R_INF=1./(1.+exp((20-svolt)/6.));
	    S_INF=1./(1.+exp((svolt+20)/5.));
	    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
	    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
	}


	#ifdef BASELINE
    D_INF=1./(1.+exp((-8-svolt)/7.5)); // original
    #endif
	#if defined (CON) || defined (CORM)
    D_INF=1./(1.+exp((-11.19-svolt)/5.26)); // control for casestudy
    #endif
    #ifdef SQT6
    D_INF=1./(1.+exp((-8.56-svolt)/6.12)); // SQT6 for casestudy
    #endif
    Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
    Bd=1.4/(1.+exp((svolt+5)/5));
    Cd=1./(1.+exp((50-svolt)/20));
    TAU_D=Ad*Bd+Cd;
    #ifdef BASELINE
    F_INF=1./(1.+exp((svolt+20)/7)); // original
    #endif
	#if defined (CON) || defined (CORM)
    F_INF=1./(1.+exp((svolt+33.90)/7.89)); // control for casestudy
    #endif
    #ifdef SQT6
    F_INF=1./(1.+exp((svolt+30.72)/8.15)); // SQT6 for casestudy
    #endif
    Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
    Bf=200./(1+exp((13-svolt)/10.));
    Cf=(180./(1+exp((svolt+30)/10)))+20;

	TAU_F=Af+Bf+Cf;
	
    /*if(svolt>0){
    	TAU_F=(Af+Bf+Cf)*2;
    }else{
    	TAU_F=Af+Bf+Cf;
    }*/
    
    F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
    Af2=600*exp(-(svolt+25)*(svolt+25)/170);
    Bf2=31/(1.+exp((25-svolt)/10));
    Cf2=16/(1.+exp((svolt+30)/10));
	
	TAU_F2=Af2+Bf2+Cf2;
	
	/*if(svolt>0){
    	TAU_F=(Af+Bf+Cf)*2;
		TAU_F2=(Af2+Bf2+Cf2)*2;
    }else{
    	TAU_F=Af+Bf+Cf;
		TAU_F2=Af2+Bf2+Cf2;
    }*/
    
    FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
    TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;
   


	//Update gates
	sm = M_INF-(M_INF-sm)*exp(-dt/TAU_M);
	sh = H_INF-(H_INF-sh)*exp(-dt/TAU_H);
	sj = J_INF-(J_INF-sj)*exp(-dt/TAU_J);
	sxr1 = Xr1_INF-(Xr1_INF-sxr1)*exp(-dt/TAU_Xr1);
	sxr2 = Xr2_INF-(Xr2_INF-sxr2)*exp(-dt/TAU_Xr2);
	sxs = Xs_INF-(Xs_INF-sxs)*exp(-dt/TAU_Xs);
	ss = S_INF-(S_INF-ss)*exp(-dt/TAU_S);
	sr = R_INF-(R_INF-sr)*exp(-dt/TAU_R);
	sd = D_INF-(D_INF-sd)*exp(-dt/TAU_D); 
	sf = F_INF-(F_INF-sf)*exp(-dt/TAU_F); 
	sf2 = F2_INF-(F2_INF-sf2)*exp(-dt/TAU_F2); 
	sfcass = FCaSS_INF-(FCaSS_INF-sfcass)*exp(-dt/TAU_FCaSS);

    //update voltage
    svolt = svolt + dt*(dvdt + dvgap_dt);
}

void TP06::setIstim(double param) {Istim = param;} // pA/pF
void TP06::setVolt(double value) {svolt = value;}
void TP06::setDt(double param) {dt = param;} // ms
void TP06::setDVgap_dt(double param) {dvgap_dt = param;}// mV/ms
double TP06::getDVgap_dt(double param) {return dvgap_dt;} // mV/ms
double TP06::getAbsDvdt() {return fabs(dvdt);} // mV/ms
double TP06::getIstim() {return Istim;} // pA/pF
double TP06::getICaL() {return ICaL;} // pA/pF
double TP06::getINa() {return INa;} // pA/pF
double TP06::getV() {return svolt;} // mV
CellType TP06::getCellType() {return ctype;}
#endif
