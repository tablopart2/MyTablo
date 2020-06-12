#ifndef BARRIEROPTION_H
#define BARRIEROPTION_H
#include "MarketParam.h"
#include "MarketParameters.h"
#include "PayoffBarrier.h"

class BarrierOption{
public:
	BarrierOption(double _refprice, signed int _expiryd,const PayoffBarrier& ThePayoff_);
	virtual ~BarrierOption();
	double Calc(MarketParam& para);
	double Calc(MarketParameters& paras);
	signed int GetExpiryd() const;
	double GetRefPrice() const;
	std::vector<double> GetResult() const;
	
protected:
	double refprice;
	signed int expiry_date;
	PayoffBarrier* ThePayoffPtr;
	std::vector<double> result;
	
};
#endif