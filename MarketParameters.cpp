#include "MarketParameters.h"
#include "Rate.h"
#include <vector>

using namespace std;

double MarketParameters::getBSVol(double t, double k) const
{
	return vol.getBSVol(t,k);
}

double MarketParameters::getTodayDivAmount(signed int n_t) const
{
	return div.getLumpsum(n_t, n_t);
}

void MarketParameters::calcLV()
{
	vol.calcLv(spot, r, q);
}

double MarketParameters::lvol(double t_axis, double s_axis) const
{
	return vol.lvol(t_axis, s_axis);
}

void MarketParameters::reset_Ivol_up()
{
	vol.reset_Ivol_up();
}

double MarketParameters::lvol_up(double t_axis, double s_axis) const
{
	return vol.lvol_up(t_axis, s_axis);
}

double MarketParameters::lvol_down(double t_axis, double s_axis) const
{
	return vol.lvol_down(t_axis, s_axis);
}

double MarketParameters::get_Lvol(int i_t, int i_k) const
{
	return vol.get_Lvol(i_t, i_k);
}

double MarketParameters::get_Lvol_hybrid(int i_t, double s_axis) const
{
	//int i_k = find_index_spot(s_axis);
	int i_k = find_index_spot2(s_axis); //SMART SEARCH 2019.12.18
	return vol.get_Lvol(i_t, i_k);
}

double MarketParameters::get_Lvol_hybrid_up(int i_t, double s_axis) const
{
	//int i_k = find_index_spot(s_axis);
	int i_k = find_index_spot2(s_axis); //SMART SEARCH 2019.12.18
	return vol.get_Lvol_up(i_t, i_k);
}

double MarketParameters::get_Lvol_hybrid_down(int i_t, double s_axis) const
{
	//int i_k = find_index_spot(s_axis);
	int i_k = find_index_spot2(s_axis); //SMART SEARCH 2019.12.18
	return vol.get_Lvol_down(i_t, i_k);
}
double MarketParameters::get_Lvol_up(int i_t, int i_k) const
{
	return vol.get_Lvol_up(i_t, i_k);

}

double MarketParameters::get_Lvol_down(int i_t, int i_k) const
{
	return vol.get_Lvol_down(i_t, i_k);
}


int MarketParameters::find_index_spot(double target) const
{
	return vol.find_index_spot(target);
}

int MarketParameters::find_index_term(double target) const
{
	return vol.find_index_term(target);
}

int MarketParameters::find_index_spot2(double target) const
{
	return vol.find_index_spot2(target);
}

void MarketParameters::print() const
{
	cout << "value date " << vdate << endl;
	cout << "spot  " << spot << endl;
	cout << "Ivol table\n";
	vol.print();
	cout << "\n";
	cout << "Lvol table\n";
	vol.printLV();
	cout << "\n";
	cout << "risk free rate \n";
	r.print();
	cout << "\n";
	cout << "dividend rate \n";
	q.print();

	
}

