#include "AsianOption.h"
#include <vector>
#include <random>
using namespace std;

double AsianOption::Calc(MarketParam & para, long numMc)
{
	//old interface don't use
	double s0 = para.get_spot();
	Rate R = para.get_rfrate();
	Rate Q = para.get_q();
	signed int vd = para.get_vdate();
	signed int expiryd = GetExpiryd();

	Vol vol = para.get_vol();
	vol.calcLv(s0, R, Q);
//	std::vector<signed int> obs
	
	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	double* mcvalues = new double[numMc];
	double* mcvalues_up = new double[numMc];
	double* mcvalues_down = new double[numMc];


	double s_tmp;
	double s_avg;
	
	signed int numObsDays = 0;

	int daydivide_ = 1;

	double* tau_p = new double[expiryd - vd + 1];
	double* r_forward_p = new double[expiryd - vd + 1];
	double* r_dc_p = new double[expiryd - vd + 1];
	double* q_forward_p = new double[expiryd - vd + 1];
	bool* is_obsDay = new bool[expiryd - vd + 1];

	for (signed int i = 0; i <= expiryd - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = R.getForward(tau_p[i]);
		r_dc_p[i] = R.getIntpRate(tau_p[i]);
		q_forward_p[i] = Q.getForward(tau_p[i]);
	}

	//initializatoin of is_obsDays
	for (int i = 0; i <= expiryd - vd; i++)
		is_obsDay[i] = false;

	//check asian obs days, obsDay vector should be ordered, check it in VBA module
	for (auto iter = obsDays.begin(); iter != obsDays.end(); iter++) {
		//cout << *iter;
		for (signed int t = vd; t <= expiryd; t++) {
			if (*iter == t) {
				is_obsDay[t-vd] = true;
				numObsDays += 1;
				break;
			}
		}
	}

	double dt = 1 / 365.0;

	for (long i = 0; i<numMc; i++)
	{
		s_tmp = s0;
		s_avg = 0.0;
	
		for (signed int t = vd; t <= expiryd; t++) {

			double short_vol = vol.lvol(tau_p[t - vd], s_tmp);
			double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
			double diff = short_vol*std::sqrt(dt);

			if (t == vd) {
				s_tmp = s_tmp;
			}else{
				s_tmp *= std::exp(drift + diff*ndist(gen));
			}

			if (is_obsDay[t - vd]) {
				s_avg += s_tmp;
			}
		}

		s_avg /= numObsDays;

		mcvalues[i] = std::exp(-r_dc_p[expiryd - vd] * tau_p[expiryd - vd])*((*ThePayoffPtr)(s_avg));
		mcvalues_up[i]= std::exp(-r_dc_p[expiryd - vd] * tau_p[expiryd - vd])*((*ThePayoffPtr)(s_avg*(1.01)));
		mcvalues_down[i] = std::exp(-r_dc_p[expiryd - vd] * tau_p[expiryd - vd])*((*ThePayoffPtr)(s_avg*(0.99)));

	}//for(i=0..)

	double npv = 0.0;
	double npv_up = 0.0;
	double npv_down = 0.0;

	for (long i = 0; i < numMc; i++) {
		npv += mcvalues[i];
		npv_up += mcvalues_up[i];
		npv_down += mcvalues_down[i];
	}

	npv /= numMc;
	npv_up /= numMc;
	npv_down /= numMc;

	result[0] = npv;
	result[1] = (npv_up - npv_down) / (s0*0.02); //delta
	result[2] = (npv_up - 2.0*npv + npv_down) / (s0*0.01) / (s0*0.01); //gamma
	result[5] = s0;

	delete[] mcvalues;
	delete[] mcvalues_up;
	delete[] mcvalues_down;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
	delete[] is_obsDay;

	return npv;
}

double AsianOption::CalcMC(MarketParameters & paras, long numMc)
{
	//paras, conti div
	signed int vd = paras.get_vdate();
	double s0 = paras.get_spot();
	paras.calcLV();

	signed int expiryd = GetExpiryd();

	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	double* mcvalues = new double[numMc];
	//double* mcvalues_up = new double[numMc];
	//double* mcvalues_down = new double[numMc];


	double s_tmp;
	double s_avg;

	signed int numObsDays = 0;

	int daydivide_ = 1;

	double* tau_p = new double[expiryd - vd + 1];
	double* r_forward_p = new double[expiryd - vd + 1];
	double* r_dc_p = new double[expiryd - vd + 1];
	double* q_forward_p = new double[expiryd - vd + 1];
	bool* is_obsDay = new bool[expiryd - vd + 1];

	for (signed int i = 0; i <= expiryd - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;
	int *idxT = new signed int[expiryd - vd + 1];
	for (int tfv = 0; tfv <= expiryd - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}
	//initializatoin of is_obsDays
	for (int i = 0; i <= expiryd - vd; i++)
		is_obsDay[i] = false;

	//check asian obs days, obsDay vector should be ordered, check it in VBA module
	for (auto iter = obsDays.begin(); iter != obsDays.end(); iter++) {
		//cout << *iter;
		for (signed int t = vd; t <= expiryd; t++) {
			if (*iter == t) {
				is_obsDay[t - vd] = true;
				numObsDays += 1;
				break;
			}
		}
	}


	for (long i = 0; i<numMc; i++)
	{
		s_tmp = s0;
		s_avg = 0.0;

		for (signed int t = vd; t <= expiryd; t++) { //start with t==vd for asian option

			double short_vol = paras.get_Lvol_hybrid(idxT[t - vd], s_tmp);
			double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
			double diff = short_vol*std::sqrt(dt);

			if (t == vd) {
				s_tmp = s_tmp; //need to save s0 to calculate average
			}
			else {
				s_tmp *= std::exp(drift + diff*ndist(gen));
			}

			if (is_obsDay[t - vd]) {
				s_avg += s_tmp;
			}
		}

		s_avg /= numObsDays;

		mcvalues[i] = std::exp(-r_dc_p[expiryd - vd] * tau_p[expiryd - vd])*((*ThePayoffPtr)(s_avg));
		//mcvalues_up[i] = std::exp(-r_dc_p[expiryd - vd] * tau_p[expiryd - vd])*((*ThePayoffPtr)(s_avg*(1.01)));
		//mcvalues_down[i] = std::exp(-r_dc_p[expiryd - vd] * tau_p[expiryd - vd])*((*ThePayoffPtr)(s_avg*(0.99)));

	}//for(i=0..)

	double npv = 0.0;
	double npv_up = 0.0;
	double npv_down = 0.0;

	for (long i = 0; i < numMc; i++) {
		npv += mcvalues[i];
		//npv_up += mcvalues_up[i];
		//npv_down += mcvalues_down[i];
	}

	npv /= numMc;
	npv_up /= numMc;
	npv_down /= numMc;

	result[0] = npv;
	//result[1] = (npv_up - npv_down) / (s0*0.02); //delta
	//result[2] = (npv_up - 2.0*npv + npv_down) / (s0*0.01) / (s0*0.01); //gamma
	//result[5] = s0;

	delete[] mcvalues;
	//delete[] mcvalues_up;
	//delete[] mcvalues_down;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
	delete[] is_obsDay;
	delete[] idxT;
	return npv;
}
