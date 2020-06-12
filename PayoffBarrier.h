#ifndef PAYOFFBARRIER_H
#define PAYOFFBARRIER_H
//#include "Payoff.h"
class PayoffBarrier {
public:
	PayoffBarrier(){};
	virtual double operator()(double spot) const=0;
	virtual void updator(double* vold, double *px, int maxi) const=0;
	virtual ~PayoffBarrier(){}
	virtual PayoffBarrier* clone() const=0;
	virtual void ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const {}

private:
	
};

#endif