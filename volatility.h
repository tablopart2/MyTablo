#ifndef VOLATILITY_H
#define VOLATILITY_H

#include <vector> 
#include "Rate.h"
using namespace std;

class Vol{
private:
	int nb_vol_term;
	int nb_vol_strike;
	mutable int init_i = 0;
	std::vector<std::vector<double> > Ivol;
	std::vector<double> vol_term;
	std::vector<double> vol_strike;
//	double **Ivol_;
//	double **Lvol_;
	std::vector<std::vector<double> > Lvol;
	std::vector<std::vector<double> > Lvol_up;
	std::vector<std::vector<double> > Lvol_down;

public:

	Vol(){}
	Vol(int nb_vol_term, int nb_vol_strike);
	Vol(vector<double>& vol_term_, vector<double>& vol_strike_, vector<vector<double> >& Ivol_);
	//real input : time, spot
	double lvol(double t_axis, double s_axis) const; 
	double lvol2(double t, double s);
	double lvol_up(double t_axis, double s_axis) const;
	double lvol_down(double t_axis, double s_axis) const; 

	void set_const_vol(double v); //initialize value with constant
	void set_vol_term(double* t, int n);
	void set_vol_strike(double* t, int n);
	void calcLv(double spot, const Rate& r, const Rate& q);
	void set_volsurface_by_term(double* row, int row_num, int nb_strike);
	void set_vol_by_point(double* arr, int nb_term, int nb_strike);
	void set_vol_by_point_vba(double* arr, int nb_term, int nb_strike);
	void set_vol_by_point_vba2(double* arr, int nb_term, int nb_strike);

	//indices input: time index, i_t, spot index
	void set_Ivol(int i_t, int i_k, double v);
	void reset_Ivol_up();
	double get_Ivol(int i_t, int i_k) const;
	double get_Lvol(int i_t, int i_k) const; 
	double get_Lvol_up(int i_t, int i_k) const;
	double get_Lvol_down(int i_t, int i_k) const;

	~Vol();
	//void setVol(char* fname);
	//double intp(double t, double s) const;
	void print() const;
	void foutLvol() const;
	void printLV() const;
	double getInpVol(double t, double k) const;
	double getBSVol(double t, double k) const;

	void Vol_up(double dv);
	int find_index_term(double t) const;
	int find_index_spot(double s) const;

	int find_index_spot2(double s) const;
};



#endif


