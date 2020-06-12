#ifndef PAYOFFCUO_H
#define PAYOFFCUO_H
#include "PayoffCall.h"
#include <string>
class PayoffCUO : public PayoffCall{
public: 
	PayoffCUO(double strike_, double upbarrier_, double rebate_);
	virtual double updator(signed int t, double* vold, double* px) const;
	//double GetStrike() const;
	double GetUpbarrier() const;
	std::string GetName() const;

	virtual ~PayoffCUO(); //dellocate Payoff*
	virtual double operator()(double spot) const;
	virtual Payoff* clone() const;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const;
	
protected:
	//double strike;
	double upbarrier;
	double rebate;
	std::string name;
};
#endif