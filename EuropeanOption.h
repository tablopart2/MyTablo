#ifndef EUROPEANOPTION_H
#define EUROPEANOPTION_H
#include "MarketParam.h"
#include "MarketParameters.h"
#include "Payoff.h"

class EuropeanOption{
public:
	EuropeanOption(double _refprice, signed int _expiryd,const Payoff& ThePayoff_);
	virtual ~EuropeanOption();
	double Calc(MarketParameters& para);
	double Calc_discrete(MarketParameters& para);
	double CalcMC(MarketParameters& paras, long numMc = 10);
	double CalcMC_discrete(MarketParameters & paras, long numMc);
	double CalcBS(MarketParameters& paras);

	void Simulation2(MarketParameters& paras, long nMC, bool db = false);
	double Simulation3(MarketParameters& paras, std::vector<double>& apath, bool db);

	signed int GetExpiryd() const;
	double GetRefPrice() const;
	std::vector<double> GetResult() const;
	void PrintResult() const;
	double BSiv(MarketParameters& paras, double p) const;


protected:
	double refprice;
	signed int expiry_date;
	Payoff* ThePayoffPtr;
	std::vector<double> result;
private:
	int getIndex(double target, double * px, int i_min, int i_max) const;
	mutable int init_i = 0;

private://**old interface, don't use it. I made it private for safety
	double Calc(MarketParam& para); //old
	double Calc2(MarketParam& para); //old
	double CalcMC(MarketParam& para, long numMc = 10);


};
#endif