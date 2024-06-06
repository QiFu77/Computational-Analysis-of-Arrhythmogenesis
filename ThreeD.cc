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
 * Last update : 07-10-2018
 */
/*
#include "SingleCell/NeonatalRatAtria.cc"
#include "SingleCell/RatSAN.cc"
#include "SingleCell/HumanVentricle.cc"
#include "SingleCell/GrandiCaMK.cc"
*/

//**********************************************************注意！在跑病理还是药理的时候要确保SingleCell中的单细胞模型处在正确的状态**********************

//#include "SingleCell/TP06.h"
#include "SingleCell/TPORd.cc"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <typeinfo>
#include "omp.h"

using namespace std;

//----- user define list -----
// multi-dimen
#define NX 102
#define NY 102
#define NZ 102 

// anisotropic
#define D1 0.1171 //0.1171 // 0.0195 // 0.1171
#define D2 D1/4.0
#define D3 D1/9.0

// spatial and temporal step(Glory. This is just not quite right)
double dx = 0.425; // 0.1; // 0.425;
double dy = 0.425; // 0.1; // 0.425;
double dz = 0.5; // 0.1176; // 0.5;
double dt = 0.005; //0.02 for human//0.1; // ms

#define OPENMP
//----- end user define list -----


/**
 * Construct geometry from geometry files
 *
 */
int*** constructGeometry()
{
	int x,y,z;// cell loop index
	// Stupid C++ can only initialize a multi-dimen array like this. If not, the array cannot be returned.
	int ***geo = new int** [NX];
	for (x = 0; x < NX; x++)
	{
		geo[x] = new int*[NY];
		for (y = 0; y < NY; y++)
		{
			geo[x][y] = new int[NZ];
			for (z = 0; z < NZ; z++)
				geo[x][y][z] = 0;
		}
	}

	printf("Construct geometry from files...\n");
	FILE *geoFile = fopen("Geometry/BensonHumanVentricle_102*102*102/heterogeneity.txt", "rt");                  //文件名已更改.
	
	
	int c;// cell type
	int count = 0;
	for (z = 0; z < NZ; z++) 
	{// be careful about the reading sequence
		for (y = 0; y < NY; y++) 
		{
			for (x = 0; x < NX; x++) 
			{ 
				fscanf(geoFile, "%d ", &c);
				if (c == 0) geo[x][y][z] = 0;	// node is outside tissue
				if (c == 1) geo[x][y][z] = 1;   // node is endocardial surface (for stimulation)
				if (c == 2) geo[x][y][z] = 2;	// node is endocardial
				if (c == 3) geo[x][y][z] = 3;	// node is midmyocardial
				if (c == 4) geo[x][y][z] = 4;	// node is epicardial

				if(x == 0 || x== NX-1 || y == 0 || y==NY-1 || z == 0 || z == NZ -1)
					if (geo[x][y][z] != 0)  cout << "Tissue exists in geometry boundary!" << endl;
				if (c != 0) count++;
			}
			fscanf(geoFile, "\n");
		}
		fscanf(geoFile, "\n");
	}
	fclose(geoFile);
	cout << "Tissue cell number = " << count << endl;
	return geo;
}


/**
 * Construct tissue according to geometry, and assign cell object for each node within tissue (effect node)
 * 
 * 
 */
Cell* *** constructTissue(int ***geo)
{
	printf("Assign tissue cell types according to geometry...\n");

	int x,y,z; // cell loop


	// Stupid C++ can only initialize a multi-dimen array like this. If not, the array cannot be returned.
	typedef Cell* Cellpt;
	
	int count = 0;

	Cellpt*** tissue = new Cellpt** [NX];
	for (x = 0; x < NX; x++)
	{
		tissue[x] = new Cellpt* [NY];
		for (y = 0; y < NY; y++)
		{
			tissue[x][y] = new Cellpt [NZ];
			for (z = 0; z < NZ; z++)
			{
				if(geo[x][y][z] > 0) // only initialize nodes within tissue
				{
					count ++;
					int celltype;
					if(geo[x][y][z] == 1) celltype = 1; // node is endocardial surface (for stimulation)
					if(geo[x][y][z] == 2) celltype = 1; // node is endocardial
					if(geo[x][y][z] == 3) celltype = 1;//2;//2; // node is midmyocardial
					if(geo[x][y][z] == 4) celltype = 1;//1;//1; // node is epicardial
					// FILE *initfile = fopen("SingleCell/NeonatalRatAtriaInitialValues.dat","r");
					// tissue[x][y][z] = new ORdHumanVentricle(celltype); // Glory for test.

					//
					
					if (geo[x][y][z] == 1 || geo[x][y][z] == 2) tissue[x][y][z] = new TPORd(ENDO);         //根据要求把它换成了TPORd单细胞模型
					if (geo[x][y][z] == 3) tissue[x][y][z] = new TPORd(MCELL);
					if (geo[x][y][z] == 4) tissue[x][y][z] = new TPORd(EPI);
					

					
				//	tissue[x][y][z] = new TP06(celltype);                                                       //类名已更改.
					tissue[x][y][z]->setDt(dt);
				}
				else 
					tissue[x][y][z] = NULL;

				// bug should be checked here.
				if (geo[x][y][z] == 0 && tissue[x][y][z] != NULL)
				{
					cout << "[x y z] = " << x << " " << y << " " << z << " type I assignment bug here." << endl;
					cout << "getV = " << tissue[x][y][z]->getV() << endl;
				}
				else if (geo[x][y][z] != 0 && tissue[x][y][z] == NULL)
					cout << "[x y z] = " << x << " " << y << " " << z << " type II assignment bug here." << endl;
			}
		}
	}


	double nn[27];
	double root2 = sqrt(2.0);
	double root3 = sqrt(3.0);
	double nix, niy, niz, imax, normalx, normaly, normalz, dislow, dis;
	int zcnt, ycnt, xcnt, isbound, nnodex, nnodey, nnodez;
	
	// script for test
	//int test1[3][3][3];
	//cout << "script for boundary test: test[-10][-1][-1] = " << test1[-10][-1][-1] << endl;

	#ifdef OPENMP
	#pragma omp parallel for private (nn, x, y, z, nix, niy, niz, imax, normalx, normaly, normalz, dislow, dis, zcnt, ycnt, xcnt, isbound, nnodex, nnodey, nnodez) schedule (static)
	#endif
	for (x = 0; x < NX; x++)
	for (y = 0; y < NY; y++)
	for (z = 0; z < NZ; z++)
	{
		if(geo[x][y][z] == 0)
		{
			// Glory. how did this code deal with out of boundary? i.e. geo[-1] 
			//	       under         on         above
			//	^     7  8  9     15 16 17     24 25 26
			//	^     4  5  6     13    14     21 22 23
			//	y     1  2  3     10 11 12     18 19 20
			//	 x>>
			/*
			nn[1]  =  geo[x-1][y-1][z-1]>0 ? 1 : 0;
			nn[2]  =  geo[x]  [y-1][z-1]>0 ? 1 : 0;
			nn[3]  =  geo[x+1][y-1][z-1]>0 ? 1 : 0;
			nn[4]  =  geo[x-1][y]  [z-1]>0 ? 1 : 0;
			nn[5]  =  geo[x]  [y]  [z-1]>0 ? 1 : 0;
			nn[6]  =  geo[x+1][y]  [z-1]>0 ? 1 : 0;
			nn[7]  =  geo[x-1][y+1][z-1]>0 ? 1 : 0;
			nn[8]  =  geo[x]  [y+1][z-1]>0 ? 1 : 0;
			nn[9]  =  geo[x+1][y+1][z-1]>0 ? 1 : 0;
			nn[10] =  geo[x-1][y-1][z]  >0 ? 1 : 0;
			nn[11] =  geo[x]  [y-1][z]  >0 ? 1 : 0;
			nn[12] =  geo[x+1][y-1][z]  >0 ? 1 : 0;
			nn[13] =  geo[x-1][y]  [z]  >0 ? 1 : 0;
			nn[14] =  geo[x+1][y]  [z]  >0 ? 1 : 0;
			nn[15] =  geo[x-1][y+1][z]  >0 ? 1 : 0;
			nn[16] =  geo[x]  [y+1][z]  >0 ? 1 : 0;
			nn[17] =  geo[x+1][y+1][z]  >0 ? 1 : 0;
			nn[18] =  geo[x-1][y-1][z+1]>0 ? 1 : 0;
			nn[19] =  geo[x]  [y-1][z+1]>0 ? 1 : 0;
			nn[20] =  geo[x+1][y-1][z+1]>0 ? 1 : 0;
			nn[21] =  geo[x-1][y]  [z+1]>0 ? 1 : 0;
			nn[22] =  geo[x]  [y]  [z+1]>0 ? 1 : 0;
			nn[23] =  geo[x+1][y]  [z+1]>0 ? 1 : 0;
			nn[24] =  geo[x-1][y+1][z+1]>0 ? 1 : 0;
			nn[25] =  geo[x]  [y+1][z+1]>0 ? 1 : 0;
			nn[26] =  geo[x+1][y+1][z+1]>0 ? 1 : 0;
			*/

			
			if(x-1 >= 0 && y-1 >= 0 && z-1 >= 0 ) nn[1]  =  geo[x-1][y-1][z-1]>0 ? 1 : 0;
			else nn[1]  =  0;

			if(y-1 >= 0 && z-1 >= 0) nn[2]  =  geo[x][y-1][z-1]>0 ? 1 : 0;
			else nn[2]  =  0;

			if(x+1 < NX && y-1 >= 0 && z-1 >= 0) nn[3]  =  geo[x+1][y-1][z-1]>0 ? 1 : 0;
			else nn[3]  =  0;

			if(x-1 >= 0 && z-1 >= 0) nn[4]  =  geo[x-1][y][z-1]>0 ? 1 : 0;
			else nn[4]  =  0;

			if(z-1 >= 0) nn[5]  =  geo[x][y][z-1]>0 ? 1 : 0;
			else nn[5]  =  0;

			if(x+1 < NX && z-1 >= 0) nn[6]  =  geo[x+1][y][z-1]>0 ? 1 : 0;
			else nn[6]  =  0;


			if(x-1 >= 0 && y+1 < NY && z-1 >= 0) nn[7]  =  geo[x-1][y+1][z-1]>0 ? 1 : 0;
			else nn[7]  =  0;


			if(y+1 < NY && z-1 >= 0) nn[8]  =  geo[x]  [y+1][z-1]>0 ? 1 : 0;
			else nn[8]  =  0;


			if(x+1 < NX && y+1 < NY && z-1 >=0) nn[9]  =  geo[x+1][y+1][z-1]>0 ? 1 : 0;
			else nn[9]  =  0;


			if(x-1 >= 0 && y-1>= 0) nn[10] =  geo[x-1][y-1][z]  >0 ?  1 : 0;
			else nn[10]  =  0;


			if(y-1 >= 0) nn[11] =  geo[x]  [y-1][z]  >0 ?  1 : 0;
			else nn[11]  =  0;


			if(x+1 < NX && y-1 >= 0) nn[12] =  geo[x+1][y-1][z]  >0 ?  1 : 0;
			else nn[12]  =  0;


			if(x-1 >= 0) nn[13] =  geo[x-1][y]  [z]  >0 ?  1 : 0;
			else nn[13]  =  0;


			if(x+1 < NX) nn[14] =  geo[x+1][y]  [z]  >0 ?  1 : 0;
			else nn[14]  =  0;


			if(x-1 >= 0 && y+1 < NY) nn[15] =  geo[x-1][y+1][z]  >0 ?  1 : 0;
			else nn[15]  =  0;


			if(y+1 < NY) nn[16] =  geo[x][y+1][z]  >0 ?  1 : 0;
			else nn[16]  =  0;


			if(x+1 < NX && y+1 < NY) nn[17] =  geo[x+1][y+1][z]  >0 ?  1 : 0;
			else nn[17]  =  0;


			if(x-1 >= 0 && y-1 >= 0 && z+1 < NZ) nn[18] =  geo[x-1][y-1][z+1]  >0 ?  1 : 0;
			else nn[18]  =  0;


			if(y-1 >= 0 && z+1 < NZ) nn[19] =  geo[x][y-1][z+1]  >0 ?  1 : 0;
			else nn[19]  =  0;


			if(x+1 < NX && y-1 >= 0 && z+1 < NZ) nn[20] =  geo[x+1][y-1][z+1]  >0 ?  1 : 0;
			else nn[20]  =  0;


			if(x-1 >= 0 && z+1 < NZ) nn[21] =  geo[x-1][y][z+1]  >0 ?  1 : 0;
			else nn[21]  =  0;


			if(z+1 < NZ) nn[22] =  geo[x]  [y]  [z+1]  >0 ?  1 : 0;
			else nn[22]  =  0;


			if(x+1 < NX && z+1 < NZ) nn[23] =  geo[x+1][y][z+1]  >0 ?  1 : 0;
			else nn[23]  =  0;


			if(x-1 >= 0 && y+1 < NY && z+1 < NZ) nn[24] =  geo[x-1][y+1][z+1]  >0 ?  1 : 0;
			else nn[24]  =  0;


			if(y+1 < NY && z+1 < NZ) nn[25] =  geo[x][y+1][z+1]  >0 ?  1 : 0;
			else nn[25]  =  0;


			if(x+1 < NX && y+1< NY && z+1 < NZ) nn[26] =  geo[x+1][y+1][z+1]  >0 ?  1 : 0;
			else nn[26]  =  0;
			
			// Glory. how did this code deal with out of boundary? i.e. geo[-1] 
			//	       under         on         above
			//	^     7  8  9     15 16 17     24 25 26
			//	^     4  5  6     13    14     21 22 23
			//	y     1  2  3     10 11 12     18 19 20
			//	 x>>
			// if (x == 29 && y == 39 && z == 55)
			// 	for (int tt = 0; tt < 27; tt++)
			// 		cout << "nn[" << tt << "] = " << nn[tt] << endl;



			// LOOP THROUGH NEAREST NEIGHBOURS TO DETERMINE HOW MANY ARE IN THE TISSUE
			isbound = 0;
			for (int count=1; count<=26; count++) if (nn[count]==1) isbound++;
			// IF AT LEAST ONE NEAREST NEIGHBOUR IS IN THE TISSUE, CALCULATE SURFACE NORMAL VECTOR FROM CURRENT FLOW BY ASSUMING TISSUE IS EQUIPOTENTIAL
			if (isbound>0) {
				// CURRENT FLOW IN X DIRECTION, WEIGHTED FOR DISTANCE
				nix =	nn[3] /root3 + nn[6] /root2 + nn[9] /root3 + 
					nn[12]/root2 + nn[14]       + nn[17]/root2 + 
					nn[20]/root3 + nn[23]/root2 + nn[26]/root3 -
					nn[1] /root3 - nn[4] /root2 - nn[7] /root3 -
					nn[10]/root2 - nn[13]       - nn[15]/root2 -
					nn[18]/root3 - nn[21]/root2 - nn[24]/root3 ;		
				// CURRENT FLOW IN Y DIRECTION, WEIGHTED FOR DISTANCE
				niy =	nn[7] /root3 + nn[8] /root2 + nn[9] /root3 +
					nn[15]/root2 + nn[16]       + nn[17]/root2 +
					nn[24]/root3 + nn[25]/root2 + nn[26]/root3 -
					nn[1] /root3 - nn[2] /root2 - nn[3] /root3 -
					nn[10]/root2 - nn[11]       - nn[12]/root2 -
					nn[18]/root3 - nn[19]/root2 - nn[20]/root3 ;		
				// CURRENT FLOW IN Z DIRECTION, WEIGHTED FOR DISTANCE
				niz =	nn[18]/root3 + nn[19]/root2 + nn[20]/root3 +
					nn[21]/root2 + nn[22]       + nn[23]/root2 +
					nn[24]/root3 + nn[25]/root2 + nn[26]/root3 -
					nn[1] /root3 - nn[2] /root2 - nn[3] /root3 -
					nn[4] /root2 - nn[5]        - nn[6] /root2 -
					nn[7] /root3 - nn[8] /root2 - nn[9] /root3 ;
				// FIND MAXIMUM OF X, Y AND Z CURRENT FLOWS FOR NORMALISATION
				imax = fabs(nix);
				if (fabs(niy) > imax) imax = fabs(niy);
				if (fabs(niz) > imax) imax = fabs(niz);
				// CALCULATE NORMALISED X, Y AND Z COMPONENTS OF THE SURFACE NORMAL VECTOR
				normalx = nix/imax;
				normaly = niy/imax;
				normalz = niz/imax;
				if (imax==0) {normalx = 0; normaly = 0; normalz = 0;}
				// FIND CLOSEST NEAREST NEIGHBOUR IN TISSUE TO SURFACE NORMAL VECTOR, THEN ASSIGN THE BOUNDARY SHELL NODE THE NEGATIVE INDEX OF THE NODE IN THE TISSUE
				// NOTE THAT NODES IN ARRAY bound[l][m][n] HAVE POSITIVE VALUE IF IN TISSUE, NEGATIVE VALUE IF IN BOUNDARY SHELL AND ZERO FOR OTHER
				dislow = 1000; // 给一个很大的初始值然后依次找最小的
				for (zcnt=-1; zcnt<=1; zcnt++)
				for (ycnt=-1; ycnt<=1; ycnt++)
				for (xcnt=-1; xcnt<=1; xcnt++) 
				{
					dis = sqrt(1.0*(xcnt-normalx)*(xcnt-normalx)+(ycnt-normaly)*(ycnt-normaly)+(zcnt-normalz)*(zcnt-normalz));
					if (x+xcnt >=0 && x+xcnt < NX && y+ycnt >= 0 && y+ycnt < NY && z+zcnt >=0 && z+zcnt < NZ )
						if (geo[x+xcnt][y+ycnt][z+zcnt] > 0 && dis<dislow) 
						{
							nnodex = xcnt;
							nnodey = ycnt;
							nnodez = zcnt;
							dislow = dis;
						}
				}
				tissue[x][y][z] = tissue[x+nnodex][y+nnodey][z+nnodez]; // 把最近的那个node的index取反存在自己里面
			} // if isbound.



			/*
			double 
			for (int a = -1; a <= 1; a++)
			for (int b = -1; b <= 1; b++)
			for (int c = -1; c <= 1; c++)
			{

				double min = 10000;
				int xbias = (x+a == NX || x+a == -1)? 0:a;
				int ybias = (y+b == NY || y+b == -1)? 0:b;
				int zbias = (z+c == NZ || z+c == -1)? 0:c;
				double distance = xbias*xbias + ybias*ybias + zbias*zbias;
				if(geo[x+xbias][y+ybias][z+zbias] > 0 && distance < min)
				{
					min = distance;
					geo[x][y][z] == -1; // boundary
					tissue[x][y][z] = tissue[x+xbias][y+ybias][z+zbias];
				}
			} */


		} 
	}
	cout << "construct tissue successful" << endl;                   //debug
	return tissue;
}


/**
 * Calculate diffusion tensor for each node within tissue (effect node)
 * 
 * a. Diffusion tensor has 9 elements but simplied to 6 due to diagonal symmetry.
 * b. More details and comments can be found in the calculation code.
 * 
 */
double**** calculateDiffTensor(int*** geo)
{
	int x,y,z,i;
	// Stupid C++ can only initialize a multi-dimen array like this. If not, the array cannot be returned.
	double ****D = new double*** [NX];
	for (x = 0; x < NX; x++)
	{
		D[x] = new double** [NY];
		for (y = 0; y < NY; y++)
		{
			D[x][y] = new double* [NZ];
			for (z = 0; z < NZ; z++)
			{
				D[x][y][z] = new double[6];
				for (i = 0; i < 6; i++)
					D[x][y][z][i] = 0;
			}
		} 
	}

	double f1,f2,f3,s1,s2,s3,n1,n2,n3;
	printf("Build in fibre orientations...\n");
	FILE *f1_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec1_x_wedge.txt", "rt"); // local system along the fibre orientation and project to global coordinate system x axis 
	FILE *f2_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec1_y_wedge.txt", "rt"); // local system along the fibre orientation and project to global coordinate system y axis
	FILE *f3_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec1_z_wedge.txt", "rt"); // local system along the fibre orientation and project to global coordinate system z axis

	FILE *s1_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec2_x_wedge.txt", "rt"); // local system "in-sheet" orientation and project to global coordinate system x axis
	FILE *s2_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec2_y_wedge.txt", "rt"); // local system "in-sheet" orientation and project to global coordinate system y axis
	FILE *s3_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec2_z_wedge.txt", "rt"); // local system "in-sheet" orientation and project to global coordinate system z axis

	FILE *n1_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec3_x_wedge.txt", "rt"); // local system sheet-vertical orientation and project to global coordinate system x axis
	FILE *n2_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec3_y_wedge.txt", "rt"); // local system sheet-vertical orientation and project to global coordinate system y axis
	FILE *n3_file = fopen ("Geometry/BensonHumanVentricle_102*102*102/vec3_z_wedge.txt", "rt"); // local system sheet-vertical orientation and project to global coordinate system z axis
	
	for (z = 0; z < NZ; z++) 
	{
		for (y = 0; y < NY; y++) 
		{
			for (x = 0; x < NX; x++) 
			{
				fscanf (f1_file, "%lf ", &f1);
				fscanf (f2_file, "%lf ", &f2);
				fscanf (f3_file, "%lf ", &f3);
				fscanf (s1_file, "%lf ", &s1);
				fscanf (s2_file, "%lf ", &s2);
				fscanf (s3_file, "%lf ", &s3);
				fscanf (n1_file, "%lf ", &n1);
				fscanf (n2_file, "%lf ", &n2);
				fscanf (n3_file, "%lf ", &n3);
				if (geo[x][y][z] > 0) // within tissue, only calculate diffusion tensor for nodes within tissue 
				{
					// Calculation steps:
					// | d11 d12 d13 |     | D[0] D[1] D[2] |        | f1*f1 f1*f2 f1*f3 |      | s1*s1 s1*s2 s1*s3 |      | n1*n1 n1*f2 n1*n3 |
					// | d21 d22 d23 | <=> | D[1] D[3] D[4] | <=> D1*| f2*f1 f2*f2 f2*f3 | + D2*| s2*s1 s2*s2 s2*s3 | + D3*| n2*n1 n2*f2 n2*n3 |
					// | d31 d32 d33 |     | D[2] D[4] D[5] |        | f3*f1 f3*f2 f3*f3 |      | s3*s1 s3*s2 s3*s3 |      | n3*n1 n3*f2 n3*n3 |
		        	D[x][y][z][0] = D1*f1*f1 + D2*s1*s1 + D3*n1*n1;
		        	D[x][y][z][1] = D1*f1*f2 + D2*s1*s2 + D3*n1*n2;
		        	D[x][y][z][2] = D1*f1*f3 + D2*s1*s3 + D3*n1*n3;
		        	D[x][y][z][3] = D1*f2*f2 + D2*s2*s2 + D3*n2*n2;
		        	D[x][y][z][4] = D1*f2*f3 + D2*s2*s3 + D3*n2*n3;
		        	D[x][y][z][5] = D1*f3*f3 + D2*s3*s3 + D3*n3*n3;

		      //   	if(x == 26 && y == 78 && (z == 15 || z == 16))
				    // {
				    // 	cout << "z = " << z << endl;
				    // 	cout << "[f1 f2 f3 s1 s2 s3 n1 n2 n3] = " << f1 << f2 << f3 << s1 << s2 << s3 << n1 << n2 << n3 << endl;
				    // 	for (int t = 0; t < 6; t++)
				    // 		cout << "D[" << t << "] = " << D[x][y][z][t] << endl;
				    // }

		        	//if(x==52 && y==52 && z==52) cout << "look at here!! D[52][52][52] = " << D[52][52][52][4] << endl;
				}
			} 
			fscanf(f1_file, "\n"); 
			fscanf(f2_file, "\n"); 
			fscanf(f3_file, "\n"); 
			fscanf(s1_file, "\n"); 
			fscanf(s2_file, "\n"); 
			fscanf(s3_file, "\n"); 
			fscanf(n1_file, "\n"); 
			fscanf(n2_file, "\n"); 
			fscanf(n3_file, "\n"); 
		} 
		fscanf(f1_file, "\n"); 
		fscanf(f2_file, "\n"); 
		fscanf(f3_file, "\n"); 
		fscanf(s1_file, "\n"); 
		fscanf(s2_file, "\n"); 
		fscanf(s3_file, "\n"); 
		fscanf(n1_file, "\n"); 
		fscanf(n2_file, "\n"); 
		fscanf(n3_file, "\n"); 
	}
	fclose (f1_file);
	fclose (f2_file);
	fclose (f3_file);
	fclose (s1_file);
	fclose (s2_file);
	fclose (s3_file);
	fclose (n1_file);
	fclose (n2_file);
	fclose (n3_file);
	return D;
}


/**
 * Calculate derivation of diffusion tensor for each node within tissue (effect node)
 *
 * 	BUILD ARRAY FOR EACH NODE, SHOWING SPATIAL RATE OF CHANGE OF DIFFUSION where 0-8 indicates a position in the array:
 *	 	 | [0] [1] [2] |     | ∂d11/∂x  ∂d12/∂x  ∂d13/∂x |     | ∂D[0]/∂x  ∂D[1]/∂x  ∂D[2]/∂x |
 *	dD = | [3] [4] [5] | <=> | ∂d21/∂y  ∂d22/∂y  ∂d23/∂y | <=> | ∂D[1]/∂y  ∂D[3]/∂y  ∂D[4]/∂y | 
 *	 	 | [6] [7] [8] |     | ∂d31/∂z  ∂d32/∂z  ∂d33/∂z |     | ∂D[2]/∂z  ∂D[4]/∂z  ∂D[5]/∂z |
 */
double**** calculateDiffDerivation(int*** geo, double**** D)
{
	printf("Calculate derivation of diffusion tensor...\n");
	
	int x,y,z; // cell loop index
	int i; // dimension loop index
	// Stupid C++ can only initialize a multi-dimen array like this. If not, the array cannot be returned.
	double ****dD = new double*** [NX];
	for (x = 0; x < NX; x++)
	{
		dD[x] = new double** [NY];
		for (y = 0; y < NY; y++)
		{
			dD[x][y] = new double* [NZ];
			for (z = 0; z < NZ; z++)
			{
				dD[x][y][z] = new double[9];
				for (i = 0; i < 9; i++)
				{
					// cout << "x = " << x << "; y = " << y << "; z = " << z << ";" <<endl;
					dD[x][y][z][i] = 0;
				}
			}
		}
	}

	// if(typeid(dD[30][38][55][1]) == typeid(double)) cout << "data type double ?!?!?!?!??!?!??!?!??!" << endl;


	// if(dD[30][38][55][1] == dD[30][38][56][1]) cout << "1start fuck!!!!!!!!" << dD[30][38][56][1] << endl;
	// if(dD[30][38][55][2] == dD[30][38][56][2]) cout << "2start fuck!!!!!!!!" << dD[30][38][56][2] << endl;
	// if(dD[30][38][55][3] == dD[30][38][56][3]) cout << "3start fuck!!!!!!!!" << dD[30][38][56][3] << endl;
	// if(dD[30][38][55][4] == dD[30][38][56][4]) cout << "4start fuck!!!!!!!!" << dD[30][38][56][4] << endl;
	// if(dD[30][38][55][5] == dD[30][38][56][5]) cout << "5start fuck!!!!!!!!" << dD[30][38][56][5] << endl;
	// if(dD[30][38][55][6] == dD[30][38][56][6]) cout << "6start fuck!!!!!!!!" << dD[30][38][56][6] << endl;
	// if(dD[30][38][55][7] == dD[30][38][56][7]) cout << "7start fuck!!!!!!!!" << dD[30][38][56][7] << endl;
	// if(dD[30][38][55][8] == dD[30][38][56][8]) cout << "8start fuck!!!!!!!!" << dD[30][38][56][8] << endl;
	

	// Neighbours' diffusion that needed in calculation. Array nbd is necessary for boundary consideration.
	// double* nbd[7]; // Glory, how to parallel this variable?
	// 	for (int i = 0; i < 7; i++)
	// 		nbd[i] = new double[6]; // Glory,还是感觉有点不对劲



	// general case without boundary
	#ifdef OPENMP
	#pragma omp parallel for private (x, y, z) schedule (static)
	#endif
	for (x = 1; x < NX-1; x++) // acording to experience, the tissue must not be on the geometry boundary
	for (y = 1; y < NY-1; y++)
	for (z = 1; z < NZ-1; z++)
	if (geo[x][y][z] != 0)
	{


		/* not bothered to consider the geometry cubox boundary
		// itself
		nbd[0]  = D[x][y][z];	

		// surface
		nbd[1]  = (x+1 == NX)? D[x][y][z] : D[x+1][y][z];
		nbd[2]  = (x-1 == -1)? D[x][y][z] : D[x-1][y][z];
		nbd[3]  = (y+1 == NY)? D[x][y][z] : D[x][y+1][z];
		nbd[4]  = (y-1 == -1)? D[x][y][z] : D[x][y-1][z];
		nbd[5]  = (z+1 == NZ)? D[x][y][z] : D[x][y][z+1];
		nbd[6]  = (z-1 == -1)? D[x][y][z] : D[x][y][z-1];

		// edge
		// nbd[7]  = (x+1 == NX || y+1 == NY)? D[x][y][z] : D[x+1][y+1][z];
		// nbd[8]  = (x+1 == NX || y-1 == -1)? D[x][y][z] : D[x+1][y-1][z];
		// nbd[9]  = (x-1 == -1 || y+1 == NY)? D[x][y][z] : D[x-1][y+1][z];
		// nbd[10] = (x-1 == -1 || y-1 == -1)? D[x][y][z] : D[x-1][y-1][z];
		// nbd[11] = (x+1 == NX || z+1 == NZ)? D[x][y][z] : D[x+1][y][z+1];
		// nbd[12] = (x+1 == NX || z-1 == -1)? D[x][y][z] : D[x+1][y][z-1];
		// nbd[13] = (x-1 == -1 || z+1 == NZ)? D[x][y][z] : D[x-1][y][z+1];
		// nbd[14] = (x-1 == -1 || z-1 == -1)? D[x][y][z] : D[x-1][y][z-1];
		// nbd[15] = (y+1 == NY || z+1 == NZ)? D[x][y][z] : D[x][y+1][z+1];
		// nbd[16] = (y+1 == NY || z-1 == -1)? D[x][y][z] : D[x][y+1][z-1];
		// nbd[17] = (y-1 == -1 || z+1 == NZ)? D[x][y][z] : D[x][y-1][z+1];
		// nbd[18] = (y-1 == -1 || z-1 == -1)? D[x][y][z] : D[x][y-1][z-1];


	    // spatial rate of change of Dxx in x direction
	    dD[x][y][z][0] = (nbd[1][0] - nbd[2][0]) / (2*dx);
	    if (geo[x+1][y][z] == 0) dD[x][y][z][0] = (nbd[0][0] - nbd[2][0]) / (dx);
	    if (geo[x-1][y][z] == 0) dD[x][y][z][0] = (nbd[1][0] - nbd[0][0]) / (dx);
	    if (geo[x+1][y][z] == 0 && geo[x-1][y][z] == 0) dD[x][y][z][0] = 0.0;

		// spatial rate of change of Dxy in x direction
	    dD[x][y][z][1] = (nbd[1][1] - nbd[2][1]) / (2*dx);
	    if (geo[x+1][y][z] == 0) dD[x][y][z][1] = (nbd[0][1] - nbd[2][1]) / (dx);
	    if (geo[x-1][y][z] == 0) dD[x][y][z][1] = (nbd[1][1] - nbd[0][1]) / (dx);
	    if (geo[x+1][y][z] == 0 && geo[x-1][y][z] == 0) dD[x][y][z][1] = 0.0;

		// spatial rate of change of Dxz in x direction
	    dD[x][y][z][2] = (nbd[1][2] - nbd[2][2]) / (2*dx);		
	    if (geo[x+1][y][z] == 0) dD[x][y][z][2] = (nbd[0][2] - nbd[2][2]) / (dx);
	    if (geo[x-1][y][z] == 0) dD[x][y][z][2] = (nbd[1][2] - nbd[0][2]) / (dx);
	    if (geo[x+1][y][z] == 0 && geo[x-1][y][z] == 0) dD[x][y][z][2] = 0.0;

		// spatial rate of change of Dyx in y direction
	    dD[x][y][z][3] = (nbd[3][1] - nbd[4][1]) / (2*dy);
	    if (geo[x][y+1][z] == 0) dD[x][y][z][3] = (nbd[0][1] - nbd[4][1]) / (dy);
	    if (geo[x][y-1][z] == 0) dD[x][y][z][3] = (nbd[3][1] - nbd[0][1]) / (dy);
	    if (geo[x][y+1][z] == 0 && geo[x][y-1][z] == 0) dD[x][y][z][3] = 0.0;

		// spatial rate of change of Dyy in y direction
	    dD[x][y][z][4] = (nbd[3][3] - nbd[4][3]) / (2*dy);		
	    if (geo[x][y+1][z] == 0) dD[x][y][z][4] = (nbd[0][3] - nbd[4][3]) / (dy);
	    if (geo[x][y-1][z] == 0) dD[x][y][z][4] = (nbd[3][3] - nbd[0][3]) / (dy);
	    if (geo[x][y+1][z] == 0 && geo[x][y-1][z] == 0) dD[x][y][z][4] = 0.0;  

		// spatial rate of change of Dyz in y direction
	    dD[x][y][z][5] = (nbd[3][4] - nbd[4][4]) / (2*dy);		
	    if (geo[x][y+1][z] == 0) dD[x][y][z][5] = (nbd[0][4] - nbd[4][4]) / (dy);
	    if (geo[x][y-1][z] == 0) dD[x][y][z][5] = (nbd[3][4] - nbd[0][4]) / (dy);
	    if (geo[x][y+1][z] == 0 && geo[x][y-1][z] == 0) dD[x][y][z][5] = 0.0;

		// spatial rate of change of Dzx in z direction
	    dD[x][y][z][6] = (nbd[5][2] - nbd[6][2]) / (2*dz);	//[z+1]-[z-1]
	    if (geo[x][y][z+1] == 0) dD[x][y][z][6] = (nbd[0][2] - nbd[6][2]) / (dz);
	    if (geo[x][y][z-1] == 0) dD[x][y][z][6] = (nbd[5][2] - nbd[0][2]) / (dz);
	    if (geo[x][y][z+1] == 0 && geo[x][y][z-1] == 0) dD[x][y][z][6] = 0.0;
	    // if(x==26 && y==78 && (z==16 || z == 15)) 
	    // {
	    // 	cout << "[x y z] = " << x << " " << y << " " << z << endl;
	    // 	cout << "geo[x][y][z+1] = " << geo[x][y][z+1] << endl;
	    // 	cout << "geo[x][y][z-1] = " << geo[x][y][z-1] << endl;
	    // 	cout << "D[x][y][z+1] = nbd[5][2] = " << nbd[5][2] << endl;
	    // 	cout << "D[x][y][z-1] = nbd[6][2] = " << nbd[6][2] << endl;
	    // 	cout << "D[x][y][z]   = nbd[0][2] = " << nbd[0][2] << endl;
	    // 	cout << "[26 78 z][6] dD = " << dD[x][y][z][6] << endl; 
	    // }

		// spatial rate of change of Dzy in z direction
	    dD[x][y][z][7] = (nbd[5][4] - nbd[6][4]) / (2*dz);		
	    if (geo[x][y][z+1] == 0) dD[x][y][z][7] = (nbd[0][4] - nbd[6][4]) / (dz);
	    if (geo[x][y][z-1] == 0) dD[x][y][z][7] = (nbd[5][4] - nbd[0][4]) / (dz);
	    if (geo[x][y][z+1] == 0 && geo[x][y][z-1] == 0) dD[x][y][z][7] = 0.0;

		// spatial rate of change of Dzz in z direction
	    dD[x][y][z][8] = (nbd[5][5] - nbd[6][5]) / (2*dz);
	    if (geo[x][y][z+1] == 0) dD[x][y][z][8] = (nbd[0][5] - nbd[6][5]) / (dz);
	    if (geo[x][y][z-1] == 0) dD[x][y][z][8] = (nbd[5][5] - nbd[0][5]) / (dz);
	    if (geo[x][y][z+1] == 0 && geo[x][y][z-1] == 0) dD[x][y][z][8] = 0.0;
		*/


 		//  	 | [0] [1] [2] |     | ∂d11/∂x  ∂d12/∂x  ∂d13/∂x |     | ∂D[0]/∂x  ∂D[1]/∂x  ∂D[2]/∂x |
 		//  dD = | [3] [4] [5] | <=> | ∂d21/∂y  ∂d22/∂y  ∂d23/∂y | <=> | ∂D[1]/∂y  ∂D[3]/∂y  ∂D[4]/∂y | 
 		//  	 | [6] [7] [8] |     | ∂d31/∂z  ∂d32/∂z  ∂d33/∂z |     | ∂D[2]/∂z  ∂D[4]/∂z  ∂D[5]/∂z |

		// spatial rate of change of Dxx in x direction
	    dD[x][y][z][0] = (D[x+1][y][z][0] - D[x-1][y][z][0]) / (2*dx);
	    if (geo[x+1][y][z] == 0) dD[x][y][z][0] = (D[x][y][z][0] - D[x-1][y][z][0]) / (dx);
	    if (geo[x-1][y][z] == 0) dD[x][y][z][0] = (D[x+1][y][z][0] - D[x][y][z][0]) / (dx);
	    if (geo[x+1][y][z] == 0 && geo[x-1][y][z] == 0) dD[x][y][z][0] = 0.0;

		// spatial rate of change of Dxy in x direction
	    dD[x][y][z][1] = (D[x+1][y][z][1] - D[x-1][y][z][1]) / (2*dx);
	    if (geo[x+1][y][z] == 0) dD[x][y][z][1] = (D[x][y][z][1] - D[x-1][y][z][1]) / (dx);
	    if (geo[x-1][y][z] == 0) dD[x][y][z][1] = (D[x+1][y][z][1] - D[x][y][z][1]) / (dx);
	    if (geo[x+1][y][z] == 0 && geo[x-1][y][z] == 0) dD[x][y][z][1] = 0.0;

		// spatial rate of change of Dxz in x direction
	    dD[x][y][z][2] = (D[x+1][y][z][2] - D[x-1][y][z][2]) / (2*dx);		
	    if (geo[x+1][y][z] == 0) dD[x][y][z][2] = (D[x][y][z][2] - D[x-1][y][z][2]) / (dx);
	    if (geo[x-1][y][z] == 0) dD[x][y][z][2] = (D[x+1][y][z][2] - D[x][y][z][2]) / (dx);
	    if (geo[x+1][y][z] == 0 && geo[x-1][y][z] == 0) dD[x][y][z][2] = 0.0;

		// spatial rate of change of Dyx in y direction
	    dD[x][y][z][3] = (D[x][y+1][z][1] - D[x][y-1][z][1]) / (2*dy);
	    if (geo[x][y+1][z] == 0) dD[x][y][z][3] = (D[x][y][z][1] - D[x][y-1][z][1]) / (dy);
	    if (geo[x][y-1][z] == 0) dD[x][y][z][3] = (D[x][y+1][z][1] - D[x][y][z][1]) / (dy);
	    if (geo[x][y+1][z] == 0 && geo[x][y-1][z] == 0) dD[x][y][z][3] = 0.0;

		// spatial rate of change of Dyy in y direction
	    dD[x][y][z][4] = (D[x][y+1][z][3] - D[x][y-1][z][3]) / (2*dy);		
	    if (geo[x][y+1][z] == 0) dD[x][y][z][4] = (D[x][y][z][3] - D[x][y-1][z][3]) / (dy);
	    if (geo[x][y-1][z] == 0) dD[x][y][z][4] = (D[x][y+1][z][3] - D[x][y][z][3]) / (dy);
	    if (geo[x][y+1][z] == 0 && geo[x][y-1][z] == 0) dD[x][y][z][4] = 0.0;  

		// spatial rate of change of Dyz in y direction
	    dD[x][y][z][5] = (D[x][y+1][z][4] - D[x][y-1][z][4]) / (2*dy);		
	    if (geo[x][y+1][z] == 0) dD[x][y][z][5] = (D[x][y][z][4] - D[x][y-1][z][4]) / (dy);
	    if (geo[x][y-1][z] == 0) dD[x][y][z][5] = (D[x][y+1][z][4] - D[x][y][z][4]) / (dy);
	    if (geo[x][y+1][z] == 0 && geo[x][y-1][z] == 0) dD[x][y][z][5] = 0.0;

		// spatial rate of change of Dzx in z direction
	    dD[x][y][z][6] = (D[x][y][z+1][2] - D[x][y][z-1][2]) / (2*dz);	//[z+1]-[z-1]
	    if (geo[x][y][z+1] == 0) dD[x][y][z][6] = (D[x][y][z][2] - D[x][y][z-1][2]) / (dz);
	    if (geo[x][y][z-1] == 0) dD[x][y][z][6] = (D[x][y][z+1][2] - D[x][y][z][2]) / (dz);
	    if (geo[x][y][z+1] == 0 && geo[x][y][z-1] == 0) dD[x][y][z][6] = 0.0;

		// spatial rate of change of Dzy in z direction
	    dD[x][y][z][7] = (D[x][y][z+1][4] - D[x][y][z-1][4]) / (2*dz);		
	    if (geo[x][y][z+1] == 0) dD[x][y][z][7] = (D[x][y][z][4] - D[x][y][z-1][4]) / (dz);
	    if (geo[x][y][z-1] == 0) dD[x][y][z][7] = (D[x][y][z+1][4] - D[x][y][z][4]) / (dz);
	    if (geo[x][y][z+1] == 0 && geo[x][y][z-1] == 0) dD[x][y][z][7] = 0.0;

		// spatial rate of change of Dzz in z direction
	    dD[x][y][z][8] = (D[x][y][z+1][5] - D[x][y][z-1][5]) / (2*dz);
	    if (geo[x][y][z+1] == 0) dD[x][y][z][8] = (D[x][y][z][5] - D[x][y][z-1][5]) / (dz);
	    if (geo[x][y][z-1] == 0) dD[x][y][z][8] = (D[x][y][z+1][5] - D[x][y][z][5]) / (dz);
	    if (geo[x][y][z+1] == 0 && geo[x][y][z-1] == 0) dD[x][y][z][8] = 0.0;
	    


	    // if(x == 30 && y == 38 && z == 55)
	    // { 
	    // 	cout << "D[z+1][5] - D[z][5] = " << D[x][y][z+1][5] - D[x][y][z][5] << endl;
	    // 	cout << "dz = " << dz << endl;
	    // 	cout << dD[x][y][z][8] << endl;
	    // }



	    // if(x == 30 && y == 38 && z == 55)
	    // {
	    // 	cout << "when z = 55: " << endl;
	    // 	cout << "geo[x+1][y][z] = " << geo[x+1][y][z] << "; geo[x-1][y][z] = " << geo[x-1][y][z] << endl;
	    // 	cout << "geo[x][y][z+1] = " << geo[x][y][z+1] << "; geo[x][y][z-1] = " << geo[x][y][z-1] << endl;
	    // 	cout << "D[x+1][y][z][0] = " << D[x+1][y][z][0] << "; D[x][y][z][0] = " << D[x][y][z][0] << endl;
		   //  cout << "D[x][y][z+1][5] = " << D[x][y][z+1][5] << "; D[x][y][z][5] = " << D[x][y][z][5] << endl;
	    // 	for (int t = 0; t < 9; t++)
	    // 		cout << "dD[" << t << "] = " << dD[x][y][z][t] << endl;
	    // }

	    // if(x == 30 && y == 38 && z == 55)
	    // { 
	    // 	cout << "when z = 55:(2rd) " << endl;
	    // 	cout << "dD[" << 8 << "] = " << dD[x][y][z][8] << endl;
	    // }


	    // if(dD[30][38][55][8] == 0 )
	    // if(x == 30 && y == 38 && z <=55 )	
	    // {
	    // 	cout << "geo[30][38][55] = " << geo[30][38][55] << endl;
	    // 	cout << "[x y z] = " << x << " " << y << " " << z << endl;
	    // 	cout << "dD[30 38 55][8] = " << dD[30][38][55][8] << endl;
	    // }



		//if(dD[30][38][55][8] == 0 )	
	 //    if(x == 30 && y == 38 && z <=60 )
	 //    {
	 //    	cout << "[x y z] = " << x << " " << y << " " << z << endl;
	 //    	cout << "dD[30 38 55][8] = " << dD[30][38][55][8] << endl;
		// }



	}// end for(XYZ & geo!=0)
	return dD;
}


/**
 * output membrane potential using .vtk format
 *
 */
void writeFile(int*** geo, Cell* *** tissue, int step, double dx, double dy, double dz)
{	
	int x,y,z;
	char* filename = new char[100];       
	sprintf(filename, "Outputs/VentricleThreeDResults%04d.vtk", step);
	FILE *datafile = fopen(filename, "wt");

	fprintf(datafile, "# vtk DataFile Version 3.0\n");
	fprintf(datafile, "vtk output\n");
	fprintf(datafile, "ASCII\n");
	fprintf(datafile, "DATASET STRUCTURED_POINTS\n");
	fprintf(datafile, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
	fprintf(datafile, "SPACING %.3f %.3f %.3f\n", dx, dy, dz);
	fprintf(datafile, "ORIGIN 0 0 0\n");
	fprintf(datafile, "POINT_DATA %d\n", NX*NY*NZ);
	fprintf(datafile, "SCALARS Voltage float 1\n");
	fprintf(datafile, "LOOKUP_TABLE default\n");


	for (z = 0; z < NZ; z++)
		for (y = 0; y < NY; y++)
		{
			for (x = 0; x < NX; x++) // Glory, the sequence should be double checked.
			{
				if (geo[x][y][z] == 0)
					fprintf (datafile, "-100 ");
				else 
					fprintf(datafile,"%.2f ", tissue[x][y][z]->getV()); // unit: mV // note that there is a space after %.2f
			}
			fprintf(datafile, "\n");
		}


	fclose(datafile);
}

int main()
{
	// parallel related stuff
	#ifdef OPENMP
	int coreNum = omp_get_num_procs();
	omp_set_num_threads(2 * coreNum);
	cout << "OpenMP is enabled in this code for parallel acceleration." << endl;
	cout << "Processor number = " << coreNum << endl;
	cout << "Threads   number = " << 2*coreNum << endl;
	#endif
	double programStartTime = omp_get_wtime();


	// user config list
	double BCL = 1000;
	double stimStrength = -20;//-60;//-20;
	double stopTime = 2*BCL; // ms
	double stimStart = 0;
	double stimDuration = 2;//ms
	FILE *datafile = fopen("Outputs/ThreeDResults.dat","w+");
	FILE *singlecell = fopen("Outputs/ThreeDSingleCellResults.dat","w+");

	// Applogies about Array:
	// Sorry i have to write them using a littel bit of confusing pointer-like type, cause stupid C++ cannot return multi-dimen array.
	// Feel free to use them like normal array. i.e. geo[x][y][z]. It's OK.

	// Rule for naming array: use strand for 1D, sheet for 2D, and tissue for 3D
	// Array 1: Tissue that consisting cell objects 
	typedef Cell* Cellpt;
	Cellpt*** tissue;

	// Array 2: Geometry that specify the tissue nodes as well as the cell types
	int*** geo;

	// Array 3: Diffusion tentsor, each node has 9 components but is simplied to 6
	//       | d11 d12 d13 |     | D[0] D[1] D[2] |
	//   D = | d21 d22 d23 | <=> | D[1] D[3] D[4] |
	//       | d31 d32 d33 |     | D[2] D[4] D[5] |
	double**** D;

	// Array 4: Derivatives of diffusion tentsor, each node has 9 components
	//	 	 | [0] [1] [2] |     | ∂d11/∂x  ∂d12/∂x  ∂d13/∂x |     | ∂D[0]/∂x  ∂D[1]/∂x  ∂D[2]/∂x |
	//  dD = | [3] [4] [5] | <=> | ∂d21/∂y  ∂d22/∂y  ∂d23/∂y | <=> | ∂D[1]/∂y  ∂D[3]/∂y  ∂D[4]/∂y | 
	//	 	 | [6] [7] [8] |     | ∂d31/∂z  ∂d32/∂z  ∂d33/∂z |     | ∂D[2]/∂z  ∂D[4]/∂z  ∂D[5]/∂z |
	double**** dD;

	// Array 5: Neighbours that involved in the diffusion
	//
	//  nb = 
	//


	// Step 1: geometry
	geo = constructGeometry();

	// Step 2: tissue (the main one in the simulation)
	tissue = constructTissue(geo);

	// Step 3: diffusion tensor
	D = calculateDiffTensor(geo);

	// Step 4: diffusion derivation
	dD = calculateDiffDerivation(geo, D);


	// Step 5: Start simulation
	double time = 0;
	double dvgap_dt = 0;
	int step = 0;


	int x, y, z; // cell loop index
	double dvdx, dvdy, dvdz, dvdx2, dvdy2, dvdz2, dvdxy, dvdxz, dvdyz; // derivatives


	printf("Environment is successfully set up! Simulation starts now.\n");
	for(time = 0.0, step = 0; time < stopTime; time += dt, step++)
	{
		if(step%50 == 0) // 50 * dt ms = 1 ms (dt = 0.02 ms)
			cout << "Simulation progress = " << 100.0*time/stopTime << "\%." << endl;

		#ifdef OPENMP
		#pragma omp parallel for private (x, y, z, dvgap_dt, dvdx, dvdy, dvdz, dvdx2, dvdy2, dvdz2, dvdxy, dvdxz, dvdyz) schedule (static)
		#endif

		for (x = 1; x < NX-1; x++) // according to experience, tissue never be on the geometry cubox boundary
		for (y = 1; y < NY-1; y++) 
		for (z = 1; z < NZ-1; z++) 
		{
			// cout << "[x y z] = [" << x << " " << y << " " << z << "] geo = " << geo[x][y][z] << endl;
			// if(x == 1 && y == 9)
				// if(tissue[1][9][20] == NULL) cout << "1 9 20 NULL" << endl;
				// else cout << "1 9 20 NOT NULL" << endl;

			if(geo[x][y][z] != 0)
			{
				// apply stimulus or not.
				tissue[x][y][z]->setIstim(0.0);
				// Apply stimulus according to user configuration (about time)
				if(time - floor(time/BCL)*BCL >= stimStart && 
			   	   time - floor(time/BCL)*BCL < stimStart + stimDuration)
				{
					if(geo[x][y][z] == 1) // endocardial outlayer cells that accept stimulus
						tissue[x][y][z]->setIstim(stimStrength);
				}
				

				// Parallel not added yet.
				// neighbours(x±1,y±1,z±1) generation, this step is necessary for the boundary consideration 
				
				// itself

				// double nb[19]; // Glory, how to parallel this variable?
				
				

				/* no need to consider geometry cubox boundary
				nb[0]  = tissue[x][y][z]->getV();				
				// x±1
				nb[1]  = (tissue[x+1][y][z] == NULL)? tissue[x][y][z]->getV() : tissue[x+1][y][z]->getV();
				nb[2]  = (tissue[x-1][y][z] == NULL)? tissue[x][y][z]->getV() : tissue[x-1][y][z]->getV();
				// y±1
				nb[3]  = (tissue[x][y+1][z] == NULL)? tissue[x][y][z]->getV() : tissue[x][y+1][z]->getV();
				nb[4]  = (tissue[x][y-1][z] == NULL)? tissue[x][y][z]->getV() : tissue[x][y-1][z]->getV();
				// z±1
				nb[5]  = (tissue[x][y][z+1] == NULL)? tissue[x][y][z]->getV() : tissue[x][y][z+1]->getV();
				nb[6]  = (tissue[x][y][z-1] == NULL)? tissue[x][y][z]->getV() : tissue[x][y][z-1]->getV();
				// x±1,y±1
				nb[7]  = (tissue[x+1][y+1][z] == NULL)? tissue[x][y][z]->getV() : tissue[x+1][y+1][z]->getV();
				nb[8]  = (tissue[x+1][y-1][z] == NULL)? tissue[x][y][z]->getV() : tissue[x+1][y-1][z]->getV();
				nb[9]  = (tissue[x-1][y+1][z] == NULL)? tissue[x][y][z]->getV() : tissue[x-1][y+1][z]->getV();
				nb[10] = (tissue[x-1][y-1][z] == NULL)? tissue[x][y][z]->getV() : tissue[x-1][y-1][z]->getV();
				// x±1,z±1
				nb[11] = (tissue[x+1][y][z+1] == NULL)? tissue[x][y][z]->getV() : tissue[x+1][y][z+1]->getV();
				nb[12] = (tissue[x+1][y][z-1] == NULL)? tissue[x][y][z]->getV() : tissue[x+1][y][z-1]->getV();
				nb[13] = (tissue[x-1][y][z+1] == NULL)? tissue[x][y][z]->getV() : tissue[x-1][y][z+1]->getV();
				nb[14] = (tissue[x-1][y][z-1] == NULL)? tissue[x][y][z]->getV() : tissue[x-1][y][z-1]->getV();
				// y±1,z±1
				nb[15] = (tissue[x][y+1][z+1] == NULL)? tissue[x][y][z]->getV() : tissue[x][y+1][z+1]->getV();
				nb[16] = (tissue[x][y+1][z-1] == NULL)? tissue[x][y][z]->getV() : tissue[x][y+1][z-1]->getV();
				nb[17] = (tissue[x][y-1][z+1] == NULL)? tissue[x][y][z]->getV() : tissue[x][y-1][z+1]->getV();
				nb[18] = (tissue[x][y-1][z-1] == NULL)? tissue[x][y][z]->getV() : tissue[x][y-1][z-1]->getV();
				*/

				// x±1
				// if(x+1 == NX) cout << "x+1 amaaaaaaazing" << endl;
				// if(x-1 == -1) cout << "x-1 amaaaaaaazing" << endl;
				// if(y+1 == NY) cout << "y+1 amaaaaaaazing" << endl;
				// if(y-1 == -1) cout << "y-1 amaaaaaaazing" << endl;
				// if(z+1 == NZ) cout << "z+1 amaaaaaaazing" << endl;
				// if(z-1 == -1) cout << "z-1 amaaaaaaazing" << endl;

				/*
				nb[1]  = tissue[x+1][y][z]->getV(); 
				nb[2]  = tissue[x-1][y][z]->getV(); 
				// y±1
				nb[3]  = tissue[x][y+1][z]->getV(); 
				nb[4]  = tissue[x][y-1][z]->getV(); 
				// z±1
				nb[5]  = tissue[x][y][z+1]->getV(); 
				nb[6]  = tissue[x][y][z-1]->getV(); 
				// x±1,y±1
				nb[7]  = tissue[x+1][y+1][z]->getV();
				nb[8]  = tissue[x+1][y-1][z]->getV();
				nb[9]  = tissue[x-1][y+1][z]->getV();
				nb[10] = tissue[x-1][y-1][z]->getV();
				// x±1,z±1
				nb[11] = tissue[x+1][y][z+1]->getV();
				nb[12] = tissue[x+1][y][z-1]->getV();
				nb[13] = tissue[x-1][y][z+1]->getV();
				nb[14] = tissue[x-1][y][z-1]->getV();
				// y±1,z±1
				nb[15] = tissue[x][y+1][z+1]->getV();
				nb[16] = tissue[x][y+1][z-1]->getV();
				nb[17] = tissue[x][y-1][z+1]->getV();
				nb[18] = tissue[x][y-1][z-1]->getV();
				*/



				// first order
				dvdx = (tissue[x+1][y][z]->getV() - tissue[x-1][y][z]->getV())/(2*dx);
				dvdy = (tissue[x][y+1][z]->getV() - tissue[x][y-1][z]->getV())/(2*dy);
				dvdz = (tissue[x][y][z+1]->getV() - tissue[x][y][z-1]->getV())/(2*dz);
				// second order
				dvdx2 = (tissue[x+1][y][z]->getV() + tissue[x-1][y][z]->getV() - 2*tissue[x][y][z]->getV())/(dx*dx);
				dvdy2 = (tissue[x][y+1][z]->getV() + tissue[x][y-1][z]->getV() - 2*tissue[x][y][z]->getV())/(dy*dy);
				dvdz2 = (tissue[x][y][z+1]->getV() + tissue[x][y][z-1]->getV() - 2*tissue[x][y][z]->getV())/(dz*dz);
				// cross deriva
				dvdxy = (tissue[x+1][y+1][z]->getV() - tissue[x+1][y-1][z]->getV()) - (tissue[x-1][y+1][z]->getV() - tissue[x-1][y-1][z]->getV())/(4*dx*dy);
				dvdxz = (tissue[x+1][y][z+1]->getV() - tissue[x+1][y][z-1]->getV()) - (tissue[x-1][y][z+1]->getV() - tissue[x-1][y][z-1]->getV())/(4*dx*dz);
				dvdyz = (tissue[x][y+1][z+1]->getV() - tissue[x][y+1][z-1]->getV()) - (tissue[x][y-1][z+1]->getV() - tissue[x][y-1][z-1]->getV())/(4*dy*dz);

				// i,j ϵ (x,y,z) or (1,2,3)
				// ΣΣ(∂dij/∂i * ∂v/∂j + dij * ∂2v/∂i∂j)
				// i = 1, j = 1,2,3  (∂d11/∂x * ∂v/∂x + d11 * ∂2v/∂x2) + (∂d12/∂x * ∂v/∂y + d12 * ∂2v/∂x∂y) + (∂d13/∂x * ∂v/∂z + d13 * ∂2v/∂x∂z)
				// i = 2, j = 1,2,3  (∂d21/∂y * ∂v/∂x + d21 * ∂2v/∂y∂x) + (∂d22/∂y * ∂v/∂y + d22 * ∂2v/∂y2) + (∂d23/∂y * ∂v/∂z + d23 * ∂2v/∂y∂z)
				// i = 3, j = 1,2,3  (∂d31/∂z * ∂v/∂x + d31 * ∂2v/∂z∂x) + (∂d32/∂z * ∂v/∂y + d32 * ∂2v/∂z∂y) + (∂d33/∂z * ∂v/∂z + d33 * ∂2v/∂x2)
				//           | d11 d12 d13 |     | D[0] D[1] D[2] |
				// And,  D = | d21 d22 d23 | <=> | D[1] D[3] D[4] |
				//           | d31 d32 d33 |     | D[2] D[4] D[5] |				
				dvgap_dt = (dD[x][y][z][0] * dvdx + D[x][y][z][0] * dvdx2) + (dD[x][y][z][1] * dvdy + D[x][y][z][1] * dvdxy) + (dD[x][y][z][2] * dvdz + D[x][y][z][2] * dvdxz) +
				           (dD[x][y][z][3] * dvdx + D[x][y][z][1] * dvdxy) + (dD[x][y][z][4] * dvdy + D[x][y][z][3] * dvdy2) + (dD[x][y][z][5] * dvdz + D[x][y][z][4] * dvdyz) +
				           (dD[x][y][z][6] * dvdx + D[x][y][z][2] * dvdxz) + (dD[x][y][z][7] * dvdy + D[x][y][z][4] * dvdyz) + (dD[x][y][z][8] * dvdz + D[x][y][z][5] * dvdz2);
				tissue[x][y][z]->setDVgap_dt(dvgap_dt);


				// update
				tissue[x][y][z]->update();
				
				// 30 38 55 error (<-100)
				// if(tissue[x][y][z]->getV() >= 70 || tissue[x][y][z]->getV() < -100)  // code should go through checking here.
				// if((x == 26 && y == 78 && z == 15) )
				
				if((x == 30 && y == 38 && z == 55))

				
				// if(x == 69 && y == 80 && z == 54)
				{
					fprintf (singlecell, "%4.10f\n", tissue[x][y][z]->getV());
					cout << "isotropic validation: " << endl;
					for (int k = 0; k < 9; k++)
						cout << "dD[" << k << "] = " << dD[x][y][z][k] << endl;
					for (int l = 0; l < 6; l++)
						cout << "D[" << l << "] = " << D[x][y][z][l] << endl;
					
					cout << "time = " << time << endl;
					cout << "stim = " << tissue[x][y][z]->getIstim() << endl;
					cout << "INa = " << tissue[x][y][z]->getINa() << endl;
					cout << "[x y z] = [" << x << " " << y << " " << z << "]" << endl;

					/*
					cout << "dvgap_dt = " << dvgap_dt << ";  tissue[x][y][z]->getDVgap_dt() = " << tissue[x][y][z]->getDVgap_dt() << endl;
					cout << "dv_dt = " << tissue[x][y][z]->getDvdt() << endl;
					cout << "delta = " << (dvgap_dt + tissue[x][y][z]->getDvdt())*dt << endl;						
					*/
					
					cout << "geo[x][y][z] = " << geo[x][y][z] << "; volt = " << tissue[x][y][z]->getV() << endl;					
					cout << "geo[x+1][y][z] = " << geo[x+1][y][z] << "; volt = " << tissue[x+1][y][z]->getV() << endl;
					cout << "geo[x-1][y][z] = " << geo[x-1][y][z] << endl; //"; volt = " << tissue[x-1][y][z]->getV() << endl;
					cout << "geo[x][y+1][z] = " << geo[x][y+1][z] << endl; //"; volt = " << tissue[x][y+1][z]->getV() << endl;
					cout << "geo[x][y-1][z] = " << geo[x][y-1][z] << "; volt = " << tissue[x][y-1][z]->getV() << endl;
					cout << "geo[x][y][z+1] = " << geo[x][y][z+1] << "; volt = " << tissue[x][y][z+1]->getV() << endl;
					cout << "geo[x][y][z-1] = " << geo[x][y][z-1] << endl; //"; volt = " << tissue[x][y][z-1]->getV() << endl;
					cout << "geo[x+1][y+1][z] = " << geo[x+1][y+1][z] << "; volt = " << tissue[x+1][y+1][z]->getV() << endl;
					cout << "geo[x+1][y-1][z] = " << geo[x+1][y-1][z] << "; volt = " << tissue[x+1][y-1][z]->getV() << endl;
					cout << "geo[x-1][y+1][z] = " << geo[x-1][y+1][z] << endl; //"; volt = " << tissue[x-1][y+1][z]->getV() << endl;					
					cout << "geo[x-1][y-1][z] = " << geo[x-1][y-1][z] << endl; //"; volt = " << tissue[x-1][y-1][z]->getV() << endl;
					cout << "geo[x+1][y][z+1] = " << geo[x+1][y][z+1] << "; volt = " << tissue[x+1][y][z+1]->getV() << endl;
					cout << "geo[x+1][y][z-1] = " << geo[x+1][y][z-1] << "; volt = " << tissue[x+1][y][z-1]->getV() << endl;
					cout << "geo[x-1][y][z+1] = " << geo[x-1][y][z+1] << endl; //"; volt = " << tissue[x-1][y][z+1]->getV() << endl;
					cout << "geo[x-1][y][z-1] = " << geo[x-1][y][z-1] << endl; //"; volt = " << tissue[x-1][y][z-1]->getV() << endl;
					cout << "geo[x][y+1][z+1] = " << geo[x][y+1][z+1] << endl; //"; volt = " << tissue[x][y+1][z+1]->getV() << endl;
					cout << "geo[x][y+1][z-1] = " << geo[x][y+1][z-1] << endl; //"; volt = " << tissue[x][y+1][z-1]->getV() << endl;
					cout << "geo[x][y-1][z+1] = " << geo[x][y-1][z+1] << "; volt = " << tissue[x][y-1][z+1]->getV() << endl;
					cout << "geo[x][y-1][z-1] = " << geo[x][y-1][z-1] << endl; //"; volt = " << tissue[x][y-1][z-1]->getV() << endl;
					cout << endl;
				}

				


			} //if(geo[x][y][z] != 0)
			// cout << "current cell loop has finished." << endl;
		}// end cell loop

		// output voltage.
		// if( floor(time/BCL) == numS1 - 1) // output final cycle only
		if (step % 50 == 0) // 50*dt = 1 ms once (dt = 0.02 ms)
		{
			writeFile(geo, tissue, step / 50, dx, dy, dz);
		}

	}// end time loop

	double programEndTime = omp_get_wtime();
	printf("All done. Running duration = %f\n",programEndTime - programStartTime);
	return 0;
}