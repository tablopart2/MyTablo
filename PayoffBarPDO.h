#pragma once
#include "PayoffBarrier.h"
#include <string>
class PayoffBarPDO : public PayoffBarrier {
public:
	PayoffBarPDO(double strike_, double downbarrier_, double rebate_);
	double GetStrike() const;
	double Getdownbarrier() const;
	std::string GetName() const;

	virtual ~PayoffBarPDO(); //dellocate Payoff*
	virtual double operator()(double spot) const;
	virtual void updator(double* vold, double* px, int maxi) const;
	virtual PayoffBarrier* clone() const;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const;

protected:
	double strike;
	double downbarrier;
	double rebate;
	std::string name;
};