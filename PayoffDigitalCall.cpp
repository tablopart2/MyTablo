#include "PayoffDigitalCall.h"
#include <algorithm>  
void PayoffDigitalCall::ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const
{
	for(int i=minnode;i<=maxnode;i++)
	{
		if(strike >=px[i] && strike <=px[i+1])
		{
			
			double tmp=std::min((strike-px[i-1])/2, (px[i+2]-strike)/2);
			px[i]=strike-tmp;
			px[i+1]=strike+tmp;
			
			dpx[i+1]=px[i+2]-px[i+1];
			dpx[i]=px[i+1]-px[i];
			dpx[i-1]=px[i]-px[i-1];
			break;
		}
	}
}
double PayoffDigitalCall::Get_strike() const
{
	return strike;
}
int PayoffDigitalCall::GetPayoffId() const
{
	return payoffId;
}
PayoffDigitalCall::PayoffDigitalCall(double strike_):strike(strike_)
{
}

double PayoffDigitalCall::operator()(double spot) const
{
	//return std::max(spot-strike, 0.0);
	return spot>strike? 1:0;
}

Payoff* PayoffDigitalCall::clone() const
{
	return new PayoffDigitalCall(*this);  //for upcasting
}
