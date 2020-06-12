#ifndef MARKETPARAM_H
#define MARKETPARAM_H
#include <vector>
#include "Rate.h"
#include <iostream>
#include "volatility.h"
using namespace std;

class MarketParam{
protected:
	signed int vdate;
	double spot;
	Vol vol;
	Rate r;
	Rate q;
public:
	MarketParam()
	{

	}

	MarketParam(signed int _vdate, double _spot, Vol _vol, Rate _r, Rate _q)
		:vdate(_vdate), spot(_spot), vol(_vol), r(_r), q(_q)
	{
	}
	double get_spot() const {return spot;}
	Vol get_vol() const {return vol;}
	const Rate& get_rfrate() const {return r;}
	const Rate& get_q() const {return q;}

	signed int get_vdate() const {return vdate;}
	
	void set_spot(double _s) {spot=_s;}
	//double get_Lvol()
	double get_Lvol(int i_t, int i_k) const;
	//void print_rate() const;
	
};
#endif