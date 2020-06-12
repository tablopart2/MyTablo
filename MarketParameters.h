#pragma once
#include <vector>
#include "Rate.h"
#include <iostream>
#include "volatility.h"
#include "Dividend.h"

using namespace std;

class MarketParameters {
protected:
	signed int vdate;
	double spot;
	Vol vol;
	Rate r;
	Rate q;
	Dividend div;
public:
	MarketParameters() {
	}

	//discrete dividend
	MarketParameters(signed int _vdate, double _spot, Vol _vol, Rate _r, Rate _q, Dividend _div)
		:vdate(_vdate), spot(_spot), vol(_vol), r(_r), q(_q), div(_div){
	}
	//continuous dividend
	MarketParameters(signed int _vdate, double _spot, Vol _vol, Rate _r, Rate _q)
		:vdate(_vdate), spot(_spot), vol(_vol), r(_r), q(_q){
	}
	
	double get_spot() const { return spot; }
	Vol get_vol() const { return vol; }
	const Rate& get_rfrate() const { return r; }
	const Rate& get_q() const { return q; }
	double getForward(double tau) const{ return r.getForward(tau); }
	double getDivForward(double tau) const { return q.getForward(tau); }
	double getIntpRate(double tau) const { return r.getIntpRate(tau); }
	double getDivIntpRate(double tau) const { return q.getIntpRate(tau); }
	double getBSVol(double t, double k) const;
	double getTodayDivAmount(signed int n_t) const;

	signed int get_vdate() const { return vdate; }

	void set_spot(double _s) { spot = _s; }
	void calcLV(); 

	double lvol(double t_axis, double s_axis) const; //search local vol from real input(time, spot)
	void reset_Ivol_up();
	double get_Lvol(int i_t, int i_k) const; //earch local vol from indices
	double get_Lvol_hybrid(int i_t, double s_axis) const;
	double get_Lvol_hybrid_up(int i_t, double s_axis) const;
	double get_Lvol_hybrid_down(int i_t, double s_axis) const;

	double lvol_up(double t_axis, double s_axis) const;
	double get_Lvol_up(int i_t, int i_k) const;

	double lvol_down(double t_axis, double s_axis) const;
	double get_Lvol_down(int i_t, int i_k) const;

	int find_index_spot(double target) const;
	int find_index_term(double target) const;

	int find_index_spot2(double target) const;
	void print() const;

};