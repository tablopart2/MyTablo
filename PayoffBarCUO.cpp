#include "PayoffBarCUO.h"
#include <algorithm>
#include <stdexcept>

void PayoffBarCUO::updator(double* vold, double* px, int maxi) const
{
	for(int i=0;i<=maxi;i++){
		if(px[i] >= upbarrier)
			vold[i]=rebate;
	}
}

PayoffBarCUO::PayoffBarCUO(double strike_, double upbarrier_, double rebate_)
	:strike(strike_), upbarrier(upbarrier_), rebate(rebate_)
{

}
double PayoffBarCUO::GetStrike() const
{
	return strike;
}
double PayoffBarCUO::GetUpbarrier() const
{
	return upbarrier;
}
std::string PayoffBarCUO::GetName() const
{
	return name;
}

PayoffBarCUO::~PayoffBarCUO()
{

}

double PayoffBarCUO::operator()(double spot) const
{
	if (spot >= upbarrier) {
		return rebate;
	}
	else {
		return std::max(spot - strike, 0.0);
	}

	throw std::logic_error("PayoffBarCUO::operator - no return value");
}
PayoffBarrier* PayoffBarCUO::clone() const
{
	return new PayoffBarCUO(*this);
}
void PayoffBarCUO::ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const
{
	for (int i = minnode; i <= maxnode; i++)
	{
		if (strike >= px[i] && strike <= px[i + 1])
		{

			double tmp = std::min((strike - px[i - 1]) / 2, (px[i + 2] - strike) / 2);
			px[i] = strike - tmp;
			px[i + 1] = strike + tmp;

			dpx[i + 1] = px[i + 2] - px[i + 1];
			dpx[i] = px[i + 1] - px[i];
			dpx[i - 1] = px[i] - px[i - 1];
			break;
		}
	}

	for (int i = minnode; i <= maxnode; i++)
	{
		if (upbarrier >= px[i] && upbarrier <= px[i + 1])
		{

			if (upbarrier - px[i] <= dpx[i] / 2.0)
				px[i] = upbarrier;
			else
				px[i + 1] = upbarrier;

			dpx[i + 1] = px[i + 2] - px[i + 1];
			dpx[i] = px[i + 1] - px[i];
			dpx[i - 1] = px[i] - px[i - 1];
			break;
		}
	}
}
