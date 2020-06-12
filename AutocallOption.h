#pragma once
#include "PayoffAutocallStd.h"
#include "MarketParam.h"
#include "MarketParameters.h"
#include "EuropeanOption.h"
#include <memory>
class AutocallOption {
public:
	AutocallOption() {}
	AutocallOption(double refprice_, signed int expiryd_, const PayoffAutocallStd& ThePayoff_, int hitflag=0);
	AutocallOption(char* csvfile);
	virtual ~AutocallOption();
	
	double Calc(MarketParameters& paras); //fdm
	double Calc_discrete(MarketParameters& paras); //fdm disc div
	double CalcMC(MarketParameters & paras, long nMC = 1000); //r,div optimized
	double CalcMC_discrete(MarketParameters & paras, long nMC = 1000); //r,div optimized
	double Simulation(MarketParameters& paras, long nMC);
	void Simulation2(MarketParameters& paras, long nMC, bool db = false);
	void Simulation2_1(MarketParameters & paras, EuropeanOption& eop, long numMC_, bool db);
	void Simulation3(MarketParameters& paras, std::vector<double>& apath, bool db = false);
	signed int GetExpiryd() const;
	double GetRefPrice() const;
	int GetHitFlag() const { return hitflag; }
	std::vector<double> GetResult() const;
	void PrintResult() const;
	void PrintSpec() const;

	//******old interface don't use
	double Calc_old(MarketParam& para); //interpolation not optimized
	double Calc(MarketParam& para);//interpolation optimized
	double CalcMC(MarketParam& para, long nMC = 1000); //r,div inerpolation not optimized
	double CalcMC_calc2(MarketParam & para, long nMC = 1000); //r,div optimized
	//** old interface

protected:
	double refprice;
	signed int expiry_date;
	PayoffAutocallStd* ThePayoffPtr;
	std::vector<double> result;
	int hitflag;
	double get_delta(double spot, double* px, double* uold, double* vold, int KIFlag,  int min, int max, int& init_i);
	double get_val(double spot, double* px, double* uold, double* vold, int KIFlag, int min, int max, int& init_i);

private:
	int getIndex(double target, double* px, int i_min, int i_max) const;
	mutable int init_i = 0;
};
