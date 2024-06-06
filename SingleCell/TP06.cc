/*
 * Ten Tusscher's human ventricle model Re-organization in C++.
 * Alternans and spiral breakup in a human ventricular tissue model.
 * 
 * Single cell stimulus protocal suggestion: -52 pA/pF * 1ms.
 * Timestep suggestion: dt = 0.2 ms.
 * Under Intellectual Property Protection.
 * 
 * 
 * Author      : Shugang Zhang <zhangshugang@hotmail.com>
 * Date        : 22-3-2019
 * Last update : 23-3-2023
 */


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


// ---NOTE: ONLY ONE OF THE FOLLOWING THREE LINES CAN BE KEPT HERE!---
//************************************四种类型的模型。
//#define WT    // 100 IKr
//#define HETE  // 50% IKr
#define DOMINATE   //25% IKr
 //#define HOMO  // 0% IKr


// ---NOTE: ONLY ONE OF THE FOLLOWING THREE LINES CAN BE KEPT HERE!---
//*************************************应用KCNQ1抗体后。
//#define ANTIBODY_30
//#define ANTIBODY_60


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

	//External concentrations
	double Ko = 5.4;
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
	double k1_ = 0.15;
	double k2_ = 0.045;
	double k3 = 0.060;
	double k4 = 0.005;//0.000015;
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
	//Parameters for Iks
	double pKNa = 0.03;
	double Gks;
	if(ctype == MCELL)  
		Gks = 0.098;
	else // EPI & ENDO
		Gks = 0.392;
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
	//Parameters for IbNa
	double GbNa = 0.00029;
	//Parameters for INaK
	double KmK = 1.0;
	double KmNa = 40.0;
	double knak = 2.724;
	//Parameters for ICaL
	double GCaL = 0.00003980;
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
	double KpCa = 0.0005;
	//Parameters for IpK;
	double GpK = 0.0146;




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


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL = GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.);
    Ito = Gto*sr*ss*(svolt-Ek);
    IKr = Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);

    #ifdef WT   // CASE STUDY
    IKr *= 1.0;
    #endif

    #ifdef HETE // CASE STUDY
    IKr *= 0.5;
    #endif

	#ifdef DOMINATE // CASE STUDY
    IKr *= 0.25;
    #endif
    #ifdef HOMO // CASE STUDY
    IKr *= 0;
    #endif

    IKs = Gks*sxs*sxs*(svolt-Eks);
    IK1 = GK1*rec_iK1*(svolt-Ek); // different from Lu
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
      IbNa  +
      ICaL  +
      IbCa  +
      INaK  +
      INaCa +
      IpCa  +
      IpK   +
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

		#ifdef ANTIBODY_30  // CASE STUDY
    Xs_INF=1./(1.+exp((-5.-svolt-10)/14.));  // 确认下这里的修改是否正确
    #endif
    #ifdef ANTIBODY_60  // CASE STUDY
    Xs_INF=1./(1.+exp((-5.-svolt-28.1)/14.));  // 确认下这里的修改是否正确
    #endif

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



    D_INF=1./(1.+exp((-8-svolt)/7.5)); // original
    Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
    Bd=1.4/(1.+exp((svolt+5)/5));
    Cd=1./(1.+exp((50-svolt)/20));
    TAU_D=Ad*Bd+Cd;
    F_INF=1./(1.+exp((svolt+20)/7)); // original
    Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
    Bf=200./(1+exp((13-svolt)/10.));
    Cf=(180./(1+exp((svolt+30)/10)))+20;
    TAU_F=Af+Bf+Cf;
    F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
    Af2=600*exp(-(svolt+25)*(svolt+25)/170);
    Bf2=31/(1.+exp((25-svolt)/10));
    Cf2=16/(1.+exp((svolt+30)/10));
    TAU_F2=Af2+Bf2+Cf2;
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