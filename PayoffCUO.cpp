#include "PayoffCUO.h"
#include <algorithm>

PayoffCUO::~PayoffCUO()
{

}
PayoffCUO::PayoffCUO(double strike_, double upbarrier_, double rebate_):PayoffCall(strike_),upbarrier(upbarrier_), rebate(rebate_)
{
}


double PayoffCUO::GetUpbarrier() const
{
	return upbarrier;
}
std::string PayoffCUO::GetName() const
{
	return name;
}

double PayoffCUO::operator()(double spot) const
{
	if(spot > upbarrier){
		return rebate;
	}else{
		return PayoffCall::operator()(spot);
	}

	throw std::logic_error("PAYOFFCUO operator()");
}
Payoff* PayoffCUO::clone() const
{
	return new PayoffCUO(*this);
}
void PayoffCUO::ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const
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

	for(int i=minnode;i<=maxnode;i++)
	{
		if(upbarrier >=px[i] && upbarrier <=px[i+1])
		{

			if(upbarrier-px[i] <=dpx[i]/2.0)
				px[i]=upbarrier;
			else
				px[i+1]=upbarrier;

			dpx[i+1]=px[i+2]-px[i+1];
			dpx[i]=px[i+1]-px[i];
			dpx[i-1]=px[i]-px[i-1];
			break;
		}
	}

}

double PayoffCUO::updator(signed int t, double* vold, double* px) const
{
	return 0;
}