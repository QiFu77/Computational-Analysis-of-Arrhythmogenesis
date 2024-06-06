/*
 * Virtual Base Class for all types of cells.
 * 
 * Under Intellectual Property Protection.
 * 
 * 
 * Author      : Shugang Zhang <zhangshugang@hotmail.com>
 * Date        : 15-09-2023
 */

#ifndef _CELL_CC_
#define _CELL_CC_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

enum CellType { LVEPI, LVENDO, EPI, ENDO, MCELL};

class Cell
{
public:
	virtual void init() {}
	virtual void init(CellType ct) {}
	virtual void setIstim(double param) {}
	virtual void setDt(double param) {}
	virtual void setDVgap_dt(double param) {}
	virtual CellType getCellType() {}
	virtual double getDVgap_dt() {}
	virtual double getV() {}
	virtual double getCai() {}
	virtual double getIstim() {}
	virtual double getINa(){}
	virtual double getINaL(){}
	virtual double getIt(){}
	virtual double getDvdt() {}
	virtual double getAbsDvdt(){} // mV/ms
	virtual double getICaL(){}
	virtual double getIKs(){}//xiugai
	virtual void update() {}
	virtual void outputAllStates(FILE *datafile) {}
	virtual void readinAllStates(FILE *datafile) {}
};
#endif
