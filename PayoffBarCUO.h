#ifndef PAYOFFBARCUO_H
#define PAYOFFBARCUO_H
#include "PayoffBarrier.h"
#include <string>
class PayoffBarCUO : public PayoffBarrier{
public: 
	PayoffBarCUO(double strike_, double upbarrier_, double rebate_);
	double GetStrike() const;
	double GetUpbarrier() const;
	std::string GetName() const;

	virtual ~PayoffBarCUO(); //dellocate Payoff*
	virtual double operator()(double spot) const;
	virtual void updator(double* vold, double* px, int maxi) const;
	virtual PayoffBarrier* clone() const;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const;
	
protected:
	double strike;
	double upbarrier;
	double rebate;
	std::string name;
};
#endif