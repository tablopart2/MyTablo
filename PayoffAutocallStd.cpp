#include "PayoffAutocallStd.h"
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;

PayoffAutocallStd::PayoffAutocallStd(int nb_autocall_, signed int* date_, double* strike_, double* coupon_, double ki_barrier_,double put_strike_,double dummy_coupon_ ,double refprice_)
	:nb_autocall(nb_autocall_), ki_barier(ki_barrier_),put_strike(put_strike_),dummy_coupon(dummy_coupon_), refprice(refprice_)
{
	
	if (date_[0] != -1 || strike_[0] != -1 || coupon_[0] != -1)
		throw std::logic_error("the first value of date, strike, coupon array should be zero");

	for (int i = 0; i <= nb_autocall; i++) {
		autocall_coupon.push_back(coupon_[i]);
		autocall_date.push_back(date_[i]);
		strike.push_back(strike_[i]);
	}
}

PayoffAutocallStd::PayoffAutocallStd(int nb_autocall_, vector<signed int> autocall_date_, vector<double> strike_, vector<double> coupon_, double ki_barrier_, double put_strike_, double dummy_coupon_, double refprice_)
	:nb_autocall(nb_autocall_),autocall_date(autocall_date_), strike(strike_), autocall_coupon(coupon_), ki_barier(ki_barrier_), put_strike(put_strike_), dummy_coupon(dummy_coupon_), refprice(refprice_)
{
}

void PayoffAutocallStd::final_updator(double * vold, double* uold, double * px,int mini, int maxi) const
{
	//knock-in
	for (int i = mini; i <= maxi; i++) {
		if (px[i] > strike.back()) {
			vold[i] = 1.0 + autocall_coupon.back();
		}
		else {
			vold[i] = 1.0 - std::max(put_strike -px[i], 0.0)/refprice;
		}
	}

	//not knocked-in
	for (int i = mini; i <= maxi; i++) {
		if (px[i] > ki_barier) {
//			uold[i] = 1.0 + autocall_coupon.back();
			uold[i] = 1.0 + dummy_coupon;
		}
		else {
			uold[i] = 1.0 - std::max(put_strike - px[i], 0.0) / refprice;
		}
	}
}

void PayoffAutocallStd::copy_v_to_u(double * v, double * u, int mini, int maxi) const
//it doesn't need to be a member of PayooffAutocallStd
{
	for (int i = mini; i <= maxi; i++) 
		u[i] = v[i];
}

void PayoffAutocallStd::updator(signed int td, double * vold, double* uold, double * px, int mini, int maxi) const
{
	auto iterX = strike.begin();
	auto iterC = autocall_coupon.begin();

	for (auto iter = autocall_date.begin(); iter != autocall_date.end(); ++iter) {
		if (*iter == td) {
			for (int i = mini; i <= maxi; i++) {
				if (px[i] > *iterX) {
					vold[i] = 1.0 + *iterC;
					uold[i] = 1.0 + *iterC;
				}
			}
			//break;
		}
		++iterX;
		++iterC;
	}
}


PayoffAutocallStd * PayoffAutocallStd::clone() const
{
	return new PayoffAutocallStd(*this);
}

void PayoffAutocallStd::ResetFDGrid(double * px, double * dpx, int minnode, int maxnode) const
{
	//for (int i = minnode; i <= maxnode; i++)
	//{
	//	if (strike >= px[i] && strike <= px[i + 1])
	//	{

	//		double tmp = std::min((strike - px[i - 1]) / 2, (px[i + 2] - strike) / 2);
	//		px[i] = strike - tmp;
	//		px[i + 1] = strike + tmp;

	//		dpx[i + 1] = px[i + 2] - px[i + 1];
	//		dpx[i] = px[i + 1] - px[i];
	//		dpx[i - 1] = px[i] - px[i - 1];
	//		break;
	//	}
	//}

	for (int i = minnode; i <= maxnode; i++)
	{
		if (ki_barier >= px[i] && ki_barier <= px[i + 1])
		{

			if (ki_barier - px[i] <= dpx[i] / 2.0)
				px[i] = ki_barier;
			else
				px[i + 1] = ki_barier;

			dpx[i + 1] = px[i + 2] - px[i + 1];
			dpx[i] = px[i + 1] - px[i];
			dpx[i - 1] = px[i] - px[i - 1];
			break;
		}
	}
}

int PayoffAutocallStd::GetNbAutocall() const
{
	return nb_autocall;
}

std::vector<signed int> PayoffAutocallStd::GetAutocall_date() const
{
	return autocall_date;
}

double PayoffAutocallStd::GetKiBarrier() const
{
	return ki_barier;
}

std::vector<double> PayoffAutocallStd::GetAutocall_strike() const
{
	return strike;
}

std::vector<double> PayoffAutocallStd::GetAutocall_coupon() const
{
	return autocall_coupon;
}

double PayoffAutocallStd::GetPutStrike() const
{
	return put_strike;
}

double PayoffAutocallStd::GetDummyCoupon() const
{
	return dummy_coupon;
}

double PayoffAutocallStd::GetRefPrice() const
{
	return refprice;
}

void PayoffAutocallStd::PrintPayoff() const
{
	cout << "payoff : ref price "  << refprice << endl;
	cout << "ki_barreir " << ki_barier << endl;
	cout << "put strike " << put_strike << endl;
	cout << "nb autocall " << nb_autocall << endl;
	cout << "dummy " << dummy_coupon << endl;
	cout << "autocall date\n";
	for (auto it = autocall_date.begin(); it != autocall_date.end(); it++)
		cout << *it << endl;
	cout << "strikes\n";
	for (auto it = strike.begin(); it != strike.end(); it++)
		cout << *it << endl;
	cout << "coupons \n";
	for (auto it = autocall_coupon.begin(); it != autocall_coupon.end(); it++)
		cout << *it << endl;

}
