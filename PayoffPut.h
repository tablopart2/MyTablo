#ifndef PAYOFFPUT_H
#define PAYOFFPUT_H
#include "Payoff.h"
#include <string>
class PayoffPut:public Payoff {
public: 
	PayoffPut(double strike_);
	double GetStrike() const;
	virtual ~PayoffPut(){} //dellocate Payoff*
	virtual double operator()(double spot) const;
	virtual Payoff* clone() const;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const;
	std::string GetName() const;
	virtual double Get_strike() const;
	virtual int GetPayoffId() const;
protected:
	double strike;
	std::string name;
};
#endif

