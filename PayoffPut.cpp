#include "PayoffPut.h"
#include <algorithm>

double PayoffPut::GetStrike() const
{
	return strike;
}
void PayoffPut::ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const
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
PayoffPut::PayoffPut(double strike_):strike(strike_)
{
	payoffId = -1;  //put -1, call +1
}

double PayoffPut::operator()(double spot) const
{
	return std::max(strike-spot, 0.0);
}

Payoff* PayoffPut::clone() const
{
	return new PayoffPut(*this);
}

std::string PayoffPut::GetName() const
{
	return name;
}

double PayoffPut::Get_strike() const
{
	return strike;
}

int PayoffPut::GetPayoffId() const
{
	return payoffId;
}
