#pragma once
#ifndef AMERICANOPTION_H
#define AMERICANOPTION_H
#include "Payoff.h"
#include "MarketParam.h"
#include "MarketParameters.h"
class AmericanOption {
public:
	AmericanOption(double _refprice, signed int _expiryd, const Payoff& ThePayoff_);
	virtual ~AmericanOption();
	double Calc(MarketParameters& paras);
	signed int GetExpiryd() const;
	double GetRefPrice() const;
	std::vector<double> GetResult() const;

	//***old interface
	double Calc(MarketParam& para);
protected:
	double refprice;
	signed int expiry_date;
	Payoff* ThePayoffPtr;
	std::vector<double> result;

};
#endif