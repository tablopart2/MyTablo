#include "PayoffBarPDO.h"

#include <algorithm>
#include <stdexcept>

void PayoffBarPDO::updator(double* vold, double* px, int maxi) const
{
	for (int i = 0; i <= maxi; i++) {
		if (px[i] <= downbarrier)
			vold[i] = rebate;
	}
}

PayoffBarPDO::PayoffBarPDO(double strike_, double downbarrier_, double rebate_)
	:strike(strike_), downbarrier(downbarrier_), rebate(rebate_)
{

}
double PayoffBarPDO::GetStrike() const
{
	return strike;
}
double PayoffBarPDO::Getdownbarrier() const
{
	return downbarrier;
}
std::string PayoffBarPDO::GetName() const
{
	return name;
}

PayoffBarPDO::~PayoffBarPDO()
{

}

double PayoffBarPDO::operator()(double spot) const
{
	if (spot <= downbarrier) {
		return rebate;
	}
	else {
		return std::max(strike - spot , 0.0);
	}

	throw std::logic_error("PayoffBarPDO::operator - no return value");
}
PayoffBarrier* PayoffBarPDO::clone() const
{
	return new PayoffBarPDO(*this);
}
void PayoffBarPDO::ResetFDGrid(double* px, double* dpx, int minnode, int maxnode) const
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
		if (downbarrier >= px[i] && downbarrier <= px[i + 1])
		{

			if (downbarrier - px[i] <= dpx[i] / 2.0)
				px[i] = downbarrier;
			else
				px[i + 1] = downbarrier;

			dpx[i + 1] = px[i + 2] - px[i + 1];
			dpx[i] = px[i + 1] - px[i];
			dpx[i - 1] = px[i] - px[i - 1];
			break;
		}
	}
}
