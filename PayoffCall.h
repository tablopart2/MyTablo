#ifndef PAYOFFCALL_H
#define PAYOFFCALL_H
#include "Payoff.h"

class PayoffCall:public Payoff{
public:
	PayoffCall(double strike_);
	virtual double operator()(double spot) const;
	virtual ~PayoffCall(){}
	virtual Payoff* clone() const;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const;
	virtual double Get_strike() const;
	virtual int GetPayoffId() const;
protected:
	double strike;
};

#endif