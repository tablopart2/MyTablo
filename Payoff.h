#ifndef PAYOFF_H
#define PAYOFF_H
#include <string>
class Payoff{
public:
	Payoff(){};
	virtual double operator()(double spot) const=0;
	virtual ~Payoff(){}
	virtual Payoff* clone() const=0;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const {}
	virtual double Get_strike() const = 0;
	virtual int GetPayoffId() const=0;
protected:
	int payoffId;
private:
	
};

#endif