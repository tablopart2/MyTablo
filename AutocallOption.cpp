#include <algorithm>
#include <random>
#include <cassert>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

#include "AutocallOption.h"
#include "EuropeanOption.h"
#include "j_fd.h"
#include "k_miscellaneous.hpp"

using namespace std;

AutocallOption::AutocallOption(double refprice_, signed int expiryd_, const PayoffAutocallStd & ThePayoff_, int hitflag_)
	:refprice(refprice_), expiry_date(expiryd_),hitflag(hitflag_)
{
	ThePayoffPtr = ThePayoff_.clone();
	result = std::vector<double>(30);
}

AutocallOption::AutocallOption(char * csvfile)
{
	ifstream infile(csvfile); // for example

	string line = "";
	string word = "";
	for (int i = 0; i < 14; i++) {
		getline(infile, line); //skip first 14 lines
	}

	getline(infile, line); //line15
	stringstream strstr(line);
	getline(strstr, word, ','); //cells(15,1)
	getline(strstr, word, ','); //cells(15,2)
	this->refprice = stod(word);
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line16
	strstr=stringstream(line);
	getline(strstr, word, ','); //cells(16,1)
	getline(strstr, word, ','); //cells(16,2)
	double spot = stod(word);
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line17
	getline(infile, line); //line18
	getline(infile, line); //line19
	strstr = stringstream(line);
	getline(strstr, word, ','); //cells(19,1)
	getline(strstr, word, ','); //cells(19,2)
	long nM = stoul(word);
	while (getline(strstr, word, ',')) {
	}

	for (int i = 0; i < 12; i++) {
		getline(infile, line); //skip first 22 lines
	}

	getline(infile, line); //line32
	strstr = stringstream(line);
	getline(strstr, word, ','); //cells(32,1)
	getline(strstr, word, ','); //cells(32,2)
	signed int vd = stoul(word);
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line33
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip first col
	getline(strstr, word, ','); //second col
	this->expiry_date = stoul(word);
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line34
	getline(infile, line); //line35
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip first col
	getline(strstr, word, ','); //second col
	int nokiflag = stoi(word);
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line36
	getline(infile, line); //line37
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip first col
	getline(strstr, word, ','); //second col
	this->hitflag = stoi(word);
	while (getline(strstr, word, ',')) {
	}


	getline(infile, line); //line38
	getline(infile, line); //line39
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip first col
	getline(strstr, word, ','); //second col
	double notional_protection = stod(word);
	double put_strike = refprice*(1.0 - notional_protection);
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line40
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip first col
	getline(strstr, word, ','); //second col
	double kilevel = stod(word);
	double kibarrier = kilevel*refprice;
	while (getline(strstr, word, ',')) {
	}

	getline(infile, line); //line41
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip first col
	getline(strstr, word, ','); //second col
	double dummycoupon = stod(word);
	while (getline(strstr, word, ',')) {
	}

	vector<signed int> vdates;
	vdates.push_back(-1);
	for (int i = 0; i < 12; i++) {
		getline(infile, line); 
		strstr = stringstream(line);
		getline(strstr, word, ','); //skip first col
		getline(strstr, word, ','); //second col
		if (!word.empty()) {
			if (stod(word) > 0)
				vdates.push_back(stoul(word));
		}
		while (getline(strstr, word, ',')) {
		}
	}

	int	nb_autocall = vdates.size()-1;

	vector<double> vstrikes;
	vstrikes.push_back(-1);

	for (int i = 0; i < 12; i++) {
		getline(infile, line);
		strstr = stringstream(line);
		getline(strstr, word, ','); //skip first col
		getline(strstr, word, ','); //second col
		if (!word.empty()) {
			if (stod(word) > 0)
				vstrikes.push_back(stod(word)*refprice);
		}
		while (getline(strstr, word, ',')) {
		}
	}

	vector<double> vcoupons;
	vcoupons.push_back(-1);
	for (int i = 0; i < 12; i++) {
		getline(infile, line);
		strstr = stringstream(line);
		getline(strstr, word, ','); //skip first col
		getline(strstr, word, ','); //second col
		if (!word.empty()) {
			if (stod(word) > 0)
				vcoupons.push_back(stod(word));
		}
		while (getline(strstr, word, ',')) {
		}
	}


	ThePayoffPtr = new PayoffAutocallStd(nb_autocall,vdates,vstrikes,vcoupons, kibarrier, put_strike,dummycoupon,refprice);
	infile.close();
	
}


AutocallOption::~AutocallOption()
{
	delete ThePayoffPtr;
}



double AutocallOption::Calc_old(MarketParam & para)
{
	double s0 = para.get_spot();
	Rate R = para.get_rfrate();
	Rate Q = para.get_q();
	signed int vd = para.get_vdate();
	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	Vol vol = para.get_vol();
	vol.calcLv(s0, R, Q);

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();
	int maxassetnodeindex = 300;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];


	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];


	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex=0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier() , 0, maxassetnodeindex);
	
	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	signed int t;
	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1);t--) {
			double tfv = (t - vd) / 365.0;
			double dt = 1 / 365.0;

			double r_forward = R.getForward(tfv);
			double q_forward = Q.getForward(tfv);

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				double short_vol = vol.lvol(tfv, px[i]);
				double short_vol_up = vol.lvol_up(tfv, px[i]);
				double short_vol_down = vol.lvol_down(tfv, px[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;
				beta[i] = (r_forward - q_forward)*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward, dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward, dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward, dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);

			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, kiindex, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, kiindex, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, kiindex, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

		}//for t

		//update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up,uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down,uold_down, px, 0, maxassetnodeindex);
		}
		if (t == vd)
			break;
	}

	double pv, pv_next, pv_up, pv_down;

	if (hitflag) { //hitted -> vold 
		pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);
	}
	else {
		pv = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, uold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, uold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, uold_down, 1, maxassetnodeindex - 1);
	}

	result.resize(30, 0.0);
	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02);
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01);
	result[3] = 0;  //vega
	
	if (vd == expiry_date) {
		result[4] = 0;
	}else{
		result[4] = pv_next - pv;  //theta
	}

	result[5] = s0;

	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	return pv;
}


double AutocallOption::Calc(MarketParam & para)
{
	//Calc2: speed improvemnet for local vol interpoation
	double s0 = para.get_spot();
	Rate R = para.get_rfrate();
	Rate Q = para.get_q();
	signed int vd = para.get_vdate();
	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	Vol vol = para.get_vol();
	vol.calcLv(s0, R, Q);

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	int maxassetnodeindex = 300;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];


	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	signed int t;
	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = R.getForward(tau_p[i]);
		r_dc_p[i] = R.getIntpRate(tau_p[i]);
		q_forward_p[i] = Q.getForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	int *idxS = new int[maxassetnodeindex + 1];
	int *idxT = new signed int[autocall_date[nb_autocall]-vd + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = vol.find_index_spot(px[i]);
	}

	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = vol.find_index_term(tfv/365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {
			//double tau = (t - vd) / 365.0; //time from vdate
			//double tfv = (t - vd) / 365.0;
			//double dt = 1 / 365.0;

			//double r_forward = R.getForward(tfv);
			//double q_forward = Q.getForward(tfv);

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				double short_vol = vol.lvol(tau_p[t-vd], px[i]);
				double short_vol_up = vol.lvol_up(tau_p[t - vd], px[i]);
				double short_vol_down = vol.lvol_down(tau_p[t - vd], px[i]);
				//double short_vol = vol.get_Lvol(idxT[t - vd], idxS[i]);
				//double short_vol_up = vol.get_Lvol_up(idxT[t - vd], idxS[i]);
				//double short_vol_down = vol.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				beta[i] = (r_forward_p[t-vd] - q_forward_p[t-vd])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t-vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t-vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t-vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);


			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
				//uold[i] = unew[i];
				//uold_up[i] = unew_up[i];
				//uold_down[i] = unew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, kiindex, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, kiindex, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, kiindex, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);


		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);

		}
		if (t == vd)
			break;

	}

	double pv, pv_next, pv_up, pv_down;

	if (hitflag) { //hitted -> vold 
		pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);
	}
	else {
		pv = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, uold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, uold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, uold_down, 1, maxassetnodeindex - 1);
	}

	result.resize(30, 0.0);
	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02);
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01);
	result[3] = 0;  //vega

	if (vd == expiry_date) {
		result[4] = 0;
	}
	else {
		result[4] = pv_next - pv;  //theta
	}

	result[5] = s0;

	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] tau_p; 
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p; 
	delete[] idxT;
	delete[] idxS;

	return pv;

}

double AutocallOption::Calc(MarketParameters & paras)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();
	
	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	int maxassetnodeindex = 400;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	signed int t;
	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0; 

	int *idxS = new int[maxassetnodeindex + 1];
	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = paras.find_index_spot(px[i]);
	}

	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				//double short_vol = paras.lvol(tau_p[t - vd], px[i]);
				//double short_vol_up = paras.lvol_up(tau_p[t - vd], px[i]);
				//double short_vol_down = paras.lvol_down(tau_p[t - vd], px[i]);
				double short_vol = paras.get_Lvol(idxT[t - vd], idxS[i]);
				double short_vol_up = paras.get_Lvol_up(idxT[t - vd], idxS[i]);
				double short_vol_down = paras.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				beta[i] = (r_forward_p[t - vd] - q_forward_p[t - vd])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);


			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, 0, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);
		}
		if (t == vd)
			break;
	}

	double pv, pv_next, pv_up, pv_down, delta1, gamma1;

	if (hitflag) { //hitted -> vold 
		pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);
		delta1 = (intp1d(s0*1.01, px, vold, 1, maxassetnodeindex - 1) - intp1d(s0*0.99, px, vold, 1, maxassetnodeindex - 1)) / (s0*1.01 - s0*0.99);
		gamma1= (intp1d(s0*1.01, px, vold, 1, maxassetnodeindex - 1) -2* intp1d(s0, px, vold, 1, maxassetnodeindex - 1) + intp1d(s0*0.99, px, vold, 1, maxassetnodeindex - 1)) / (s0*s0*0.01*0.01);
	}
	else {
		pv = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, uold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, uold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, uold_down, 1, maxassetnodeindex - 1);
		delta1 = (intp1d(s0*1.01, px, uold, 1, maxassetnodeindex - 1) - intp1d(s0*0.99, px, uold, 1, maxassetnodeindex - 1)) / (s0*1.01 - s0*0.99);
		gamma1 = (intp1d(s0*1.01, px, uold, 1, maxassetnodeindex - 1) - 2 * intp1d(s0, px, uold, 1, maxassetnodeindex - 1) + intp1d(s0*0.99, px, uold, 1, maxassetnodeindex - 1)) / (s0*s0*0.01*0.01);
	}

	result.resize(30, 0.0);
	for (auto iter = result.begin(); iter != result.end(); iter++)
		*iter = 0.0;

	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02);
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01);
	
	if (vd == expiry_date) {
		result[4] = 0;
	}
	else {
		result[4] = pv_next - pv;  //theta
	}

	result[6] = delta1;
	result[7] = gamma1;

	//pure delta


	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;

	delete[] idxT;
	delete[] idxS;

	return pv;
}

double AutocallOption::Calc_discrete(MarketParameters & paras)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();

	int nb_autocall = ThePayoffPtr->GetNbAutocall();
	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	int maxassetnodeindex = 400;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	signed int t;
	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	double dt = 1 / 365.0;

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		//q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}
	for (signed int t = autocall_date[nb_autocall]; t >= vd; t--) {
		//CAUTION: t instead of tau
		q_forward_p[t - vd] = paras.getTodayDivAmount(t) / s0 / dt;
	}

	int *idxS = new int[maxassetnodeindex + 1];
	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = paras.find_index_spot(px[i]);
	}

	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				double short_vol = paras.get_Lvol(idxT[t - vd], idxS[i]);
				double short_vol_up = paras.get_Lvol_up(idxT[t - vd], idxS[i]);
				double short_vol_down = paras.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				//CAUTION2: check calculation day t-1 is exdate or not(t is thd 1 day after cal day) 
				//check when t=vd+1, t-vd-1==0
				beta[i] = (r_forward_p[t - vd] - q_forward_p[t - vd-1])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);


			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, 0, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);
		}
		if (t == vd)
			break;
	}

	double pv, pv_next, pv_up, pv_down, delta1, gamma1;

	if (hitflag) { //hitted -> vold 
		pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);
		delta1 = (intp1d(s0*1.01, px, vold, 1, maxassetnodeindex - 1) - intp1d(s0*0.99, px, vold, 1, maxassetnodeindex - 1)) / (s0*1.01 - s0*0.99);
		gamma1 = (intp1d(s0*1.01, px, vold, 1, maxassetnodeindex - 1) - 2 * intp1d(s0, px, vold, 1, maxassetnodeindex - 1) + intp1d(s0*0.99, px, vold, 1, maxassetnodeindex - 1)) / (s0*s0*0.01*0.01);
	}
	else {
		pv = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, uold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, uold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, uold_down, 1, maxassetnodeindex - 1);
		delta1 = (intp1d(s0*1.01, px, uold, 1, maxassetnodeindex - 1) - intp1d(s0*0.99, px, uold, 1, maxassetnodeindex - 1)) / (s0*1.01 - s0*0.99);
		gamma1 = (intp1d(s0*1.01, px, uold, 1, maxassetnodeindex - 1) - 2 * intp1d(s0, px, uold, 1, maxassetnodeindex - 1) + intp1d(s0*0.99, px, uold, 1, maxassetnodeindex - 1)) / (s0*s0*0.01*0.01);
	}

	result.resize(30, 0.0);
	for (auto iter = result.begin(); iter != result.end(); iter++)
		*iter = 0.0;

	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02);
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01);

	if (vd == expiry_date) {
		result[4] = 0;
	}
	else {
		result[4] = pv_next - pv;  //theta
	}

	result[6] = delta1;
	result[7] = gamma1;

	//pure delta


	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;

	delete[] idxT;
	delete[] idxS;

	return pv;
}
double AutocallOption::Simulation(MarketParameters & paras, long numMC_)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();

	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	int maxassetnodeindex = 400;

	double **vgrid;
	double **ugrid;
	int nb_date = autocall_date[nb_autocall] - vd+1;
	vgrid = new double*[nb_date];
	ugrid = new double*[nb_date];
	for (int tau = 0; tau < nb_date; tau++) {
		vgrid[tau] = new double[maxassetnodeindex + 1];
		ugrid[tau] = new double[maxassetnodeindex + 1];
	}

	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	//final u,v
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vgrid[nb_date - 1][i] = vold[i];
		ugrid[nb_date - 1][i] = uold[i];
	}
	

	signed int t;
	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	int *idxS = new int[maxassetnodeindex + 1];
	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = paras.find_index_spot(px[i]);
	}

	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				//double short_vol = paras.lvol(tau_p[t - vd], px[i]);
				//double short_vol_up = paras.lvol_up(tau_p[t - vd], px[i]);
				//double short_vol_down = paras.lvol_down(tau_p[t - vd], px[i]);
				double short_vol = paras.get_Lvol(idxT[t - vd], idxS[i]);
				double short_vol_up = paras.get_Lvol_up(idxT[t - vd], idxS[i]);
				double short_vol_down = paras.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				beta[i] = (r_forward_p[t - vd] - q_forward_p[t - vd])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);


			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, 0, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

			for (int i = 0; i <= maxassetnodeindex; i++) {
				vgrid[t-vd - 1][i] = vold[i];
				ugrid[t-vd - 1][i] = uold[i];
			}
		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);
			
			//save vold, uold after updating
			//for (int i = 0; i <= maxassetnodeindex; i++) {
			//	vgrid[t - vd - 1][i] = vold[i];
			//	ugrid[t - vd - 1][i] = uold[i];
			//}

		}
		if (t == vd)
			break;
	}

	double pv, pv_next, pv_up, pv_down;

	if (hitflag) { //hitted -> vold 
		pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);
	}
	else {
		pv = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
		pv_next = intp1d(s0, px, uold_next, 1, maxassetnodeindex - 1);
		pv_up = intp1d(s0*1.01, px, uold_up, 1, maxassetnodeindex - 1);
		pv_down = intp1d(s0*0.99, px, uold_down, 1, maxassetnodeindex - 1);
	}

	result.resize(30, 0.0);
	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02);
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01);
	result[3] = 0;  //vega

	if (vd == expiry_date) {
		result[4] = 0;
	}
	else {
		result[4] = pv_next - pv;  //theta
	}

	result[5] = s0;

	//MC simulation
	double kibarrier = ThePayoffPtr->GetKiBarrier();
	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();
	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();
	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	double* mcvalues = new double[numMC_];
	vector<vector<double> > paths;
	vector<vector<double> > deltas;
	vector<vector<double> > vals;
	vector<vector<double> > cashes;
	vector <vector<double> > npvs;

	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	for (long i = 0; i<numMC_; i++)
	{
		vector<double> apath;
		vector<double> adelta;
		vector<double> aval;
		vector<double> acash;
		vector<double> anpv;

		int init_i = 0;
		int tau = 0;
		s_tmp = s0;
		tmpKIFlag = hitflag;
		
		apath.push_back(s0);
		//init_i = get_spot_index(s0, px, 0, maxassetnodeindex, init_i);
		adelta.push_back(get_delta(s0, px, ugrid[0], vgrid[0], tmpKIFlag, 0, maxassetnodeindex, init_i));
		aval.push_back(get_val(s0,px,ugrid[0], vgrid[0], tmpKIFlag, 0, maxassetnodeindex, init_i));
		acash.push_back(aval.back() - adelta.back()*s0);
		anpv.push_back(acash.back()-aval.back()+adelta.back()*s0); //should be zero
		
		assert(round(anpv[0]*100000.0)/100000.0 == 0.0);
		assert(aval.front() == pv);


		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				//double short_vol = paras.lvol(tau_p[t - vd], s_tmp);
				double short_vol = paras.get_Lvol_hybrid(idxT[t - vd], s_tmp);

				double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
				if (s_tmp<kibarrier)
					tmpKIFlag = 1;

				double df = std::exp(-r_dc_p[t - vd] * tau_p[t - vd]);

				apath.push_back(s_tmp);
				adelta.push_back(get_delta(s_tmp, px, ugrid[t-vd], vgrid[t-vd], tmpKIFlag, 0, maxassetnodeindex, init_i));
				aval.push_back(df*(get_val(s_tmp, px, ugrid[t - vd], vgrid[t - vd], tmpKIFlag, 0, maxassetnodeindex, init_i)));
				acash.push_back(df*(acash.back()*std::exp(r_forward_p[t - vd] * dt) - s_tmp*(adelta.back() - *(adelta.rbegin() + 1))));

				anpv.push_back(acash.back() - aval.back() + df*s_tmp*adelta.back());
			}

			if (s_tmp >= autocall_strike[k]) { //check autocallability
				double df = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd]);
				mcvalues[i] = df*(1.0 + autocall_coupon[k]);
				aval.back() = mcvalues[i]; 
				anpv.back()=acash.back() - aval.back() + df*s_tmp*adelta.back();
				break; //k loop
			}

			//we are here because it hasn't been autocalled
			if (k == nb_autocall) {
				double df = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd]);
				if (s_tmp >= autocall_strike[k]) {
					
					mcvalues[i] = df*(1.0 + autocall_coupon[k]);
					aval.back() = mcvalues[i];
					anpv.back() = acash.back() - aval.back() + df*s_tmp*adelta.back();
				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues[i] =df*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
						aval.back() = mcvalues[i];
						anpv.back() = acash.back() - aval.back() + df*s_tmp*adelta.back();

					}
					else if (tmpKIFlag == 0) {
						mcvalues[i] = df*(1.0 + dummy_coupon);
						aval.back() = mcvalues[i];
						anpv.back() = acash.back() - aval.back() + df*s_tmp*adelta.back();
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					aval.back() = mcvalues[i];
					anpv.back() = acash.back() - aval.back() + s_tmp*adelta.back();
				}

			} //if k

		}//for k

		paths.push_back(apath);
		deltas.push_back(adelta);
		vals.push_back(aval);
		cashes.push_back(acash);
		npvs.push_back(anpv);

	}//for(i=0..)



	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;

	delete[] px;  //get_delta에서 쓰임
	delete[] dpx;
	delete[] idxT; //get_delta
	delete[] idxS;

	vector<double> fnpv;
	for (auto i = npvs.begin(); i != npvs.end(); i++) {
		cout << (*i).back() << endl;
		fnpv.push_back((*i).back()); //final npvs
	}

	ofstream ofs;
	ofs.open("out.txt");
	ofs.clear();
	ofs.close();
	return 0;
}

void AutocallOption::Simulation2(MarketParameters & paras, long numMC_, bool db)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();

	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	paras.calcLV();


	//if (db) {//save local vol in file
	//	ofstream flvol("Lvol.csv");
	//	//fuold << "PX, uold" << endl;
	//	for (int i = 0; i <10; i++) {
	//		for (int j = 0; j < 17; j++) {
	//			flvol << paras.get_Lvol(i, j) << ",";
	//		}
	//		flvol << endl;
	//	}
	//	flvol.close();
	//}

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	vector<vector<double> > paths;

	int hitFlag = this->hitflag;
	double refprice = this->refprice;

	assert(refprice == paras.get_spot());
	assert(paras.get_vdate() == vd);

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	vector<double> mcvalues;
	vector<double> pvs;
	vector<vector<double> > PLs;

	int maxassetnodeindex = 400;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	/*save vold, uold in vector*/
	vector<vector<double> > vgrid;
	vector<vector<double> > ugrid;

	vgrid.push_back(vector<double>(vold, vold + (maxassetnodeindex + 1)));
	ugrid.push_back(vector<double>(uold, uold + (maxassetnodeindex + 1)));

	signed int t;
	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	int *idxS = new int[maxassetnodeindex + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = paras.find_index_spot(px[i]);
	}

	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				double short_vol = paras.get_Lvol(idxT[t - vd], idxS[i]);
				double short_vol_up = paras.get_Lvol_up(idxT[t - vd], idxS[i]);
				double short_vol_down = paras.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				beta[i] = (r_forward_p[t - vd] - q_forward_p[t - vd])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);

			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, 0, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

			vgrid.push_back(vector<double>(vold, vold + (maxassetnodeindex + 1)));
			ugrid.push_back(vector<double>(uold, uold + (maxassetnodeindex + 1)));

		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);
			//업데이트후의 값으로 저장하기
			vgrid.back() = (vector<double>(vold, vold + (maxassetnodeindex + 1)));
			ugrid.back() = (vector<double>(uold, uold + (maxassetnodeindex + 1)));
		}
		if (t == vd)
			break;
	}

	double pv_fd;
	pv_fd = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
	//doule check
	auto it_vgrid = vgrid.rbegin();
	auto it_ugrid = ugrid.rbegin();

	pv_fd = intp1d(s0, px, *it_ugrid, 0, maxassetnodeindex, 0);

	cout << "npv_fd " << pv_fd << endl;

	string fn = getFnameTimeStartingWith(string("ts"));
	ofstream fout_ts(fn.c_str());
	fout_ts << "tau,s_tmp,cash,delta,pv,r_for,q,PL" << endl;

	for (long i = 0; i<numMC_; i++)
	{
		vector<double> path;
		vector<double> aPL;
		path.push_back(paras.get_spot());

		auto riter_vgrid = vgrid.rbegin();
		auto riter_ugrid = ugrid.rbegin();

		s_tmp = s0;
		tmpKIFlag = hitFlag;

		int init_i = 0;
		double cash = 0.0;
		double pv = 0.0;
		double delta = 0.0;
		double delta_new = 0.0;
		double PL = 0.0;

		int spot_idx = getIndex(s_tmp, px, 0, maxassetnodeindex);

		if (spot_idx == 0) {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);

			}
			else {
				delta = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}
		}
		else if (spot_idx == maxassetnodeindex) {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
			}
			else {
				delta = ((*riter_ugrid)[spot_idx] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}

		}
		else {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
			}
			else {
				delta = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}
		}

		cash += pv - s_tmp*delta;
		PL = cash - pv + s_tmp*delta;
		aPL.push_back(PL);
		if (db)
			fout_ts << 0 << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				double short_vol = paras.get_Lvol_hybrid(idxT[t - vd], s_tmp);
				double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				for (long t2 = 1; t2 <= daydivide_; t2++) {
					s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
					path.push_back(s_tmp);

					if (s_tmp<kibarrier)
						tmpKIFlag = 1;
				}

				spot_idx = getIndex(s_tmp, px, 0, maxassetnodeindex);

				riter_vgrid++;
				riter_ugrid++;

				if (spot_idx == 0) {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}
				}
				else if (spot_idx == maxassetnodeindex) {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}

				}
				else {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}
				}

				cash *= std::exp(r_forward_p[t - vd] * dt);
				cash -= s_tmp*(delta_new - delta);
				cash += s_tmp*delta*(std::exp(q_forward_p[t - vd] * dt) - 1.0);

				delta = delta_new;
				PL = cash - pv + s_tmp*delta;
				aPL.push_back(PL);
				if (db)
					fout_ts << t - vd << "," << s_tmp << "," << cash << "," << delta << "," << pv << "," << r_forward_p[t - vd]  << "," << q_forward_p[t - vd] <<","<< PL << endl;
			}

			//autocall payoff 
			double df = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd]);

			if (s_tmp >= autocall_strike[k]) { //check autocallability
											   //mcvalues.push_back(std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]));
				mcvalues.push_back(1.0 + autocall_coupon[k]);

				//pv = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(cash + s_tmp*delta)-mcvalues.back();
				pv = mcvalues.back();

				pvs.push_back(pv);
				PL = cash - pv + s_tmp*delta;
				aPL.back() = PL;
				if (db)
					fout_ts << tmpKIFlag << "_" << k << "-th Autocalled i=" << i <<"," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;
				break; //k loop
			}

			if (k == nb_autocall) {//we are here because it hasn't been autocalled
				if (s_tmp >= autocall_strike[k]) {
					//mcvalues.push_back(std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]));
					mcvalues.push_back(1.0 + autocall_coupon[k]);

					pv = mcvalues.back();
					pvs.push_back(pv);
					PL = cash - pv + s_tmp*delta;
					aPL.back() = PL;
					if (db)
						fout_ts << tmpKIFlag << "_" << k << "-th Autocalled i=" << i << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues.push_back(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
						pv = mcvalues.back();
						pvs.push_back(pv);
						PL = cash - pv + s_tmp*delta;
						aPL.back() = PL;
						if (db)
							fout_ts << tmpKIFlag << "_" << k << "-th Autocalled i=" << i << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

					}
					else if (tmpKIFlag == 0) {
						mcvalues.push_back(1.0 + dummy_coupon);
						pv = mcvalues.back();
						pvs.push_back(pv);
						PL = cash - pv + s_tmp*delta;
						aPL.back() = PL;
						if (db)
							fout_ts << tmpKIFlag << "_" << k << "-th Autocalled i=" << i << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues.push_back(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					pv = mcvalues.back();
					pvs.push_back(pv);
					PL = cash - pv + s_tmp*delta;
					aPL.back() = PL;
					if (db)
						fout_ts << tmpKIFlag << "_" << k << "-th Autocalled i=" << i << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;
				}

			} //if k

		}//for k

		paths.push_back(path);
		PLs.push_back(aPL);

	}//for(i=0..)

	fout_ts.close();

	string fn2 = getFnameTimeStartingWith(string("PL"));
	ofstream fout(fn2.c_str());
	auto it = PLs.begin();
	fout << "numMC,PL,callday" << endl;
	for (auto iter = PLs.begin(); iter != PLs.end(); iter++)
		fout << iter - it << "," << (*iter).back() <<","<<(*iter).size() << endl;

	fout.close();

	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] idxT;
	delete[] idxS;
	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
}

void AutocallOption::Simulation2_1(MarketParameters & paras, EuropeanOption& eop,  long numMC_, bool db)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();

	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	vector<vector<double> > paths;

	int hitFlag = this->hitflag;
	double refprice = this->refprice;

	assert(refprice == paras.get_spot());
	assert(paras.get_vdate() == vd);

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	vector<double> mcvalues;
	vector<double> pvs;
	vector<vector<double> > PLs;
	
	vector<double> pvs_final_euro;

	int maxassetnodeindex = 400;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}

	/*save vold, uold in vector*/
	vector<vector<double> > vgrid;
	vector<vector<double> > ugrid;

	vgrid.push_back(vector<double>(vold, vold + (maxassetnodeindex + 1)));
	ugrid.push_back(vector<double>(uold, uold + (maxassetnodeindex + 1)));

	signed int t;
	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	int *idxS = new int[maxassetnodeindex + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = paras.find_index_spot(px[i]);
	}

	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				double short_vol = paras.get_Lvol(idxT[t - vd], idxS[i]);
				double short_vol_up = paras.get_Lvol_up(idxT[t - vd], idxS[i]);
				double short_vol_down = paras.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				beta[i] = (r_forward_p[t - vd] - q_forward_p[t - vd])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);

			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, 0, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

			vgrid.push_back(vector<double>(vold, vold + (maxassetnodeindex + 1)));
			ugrid.push_back(vector<double>(uold, uold + (maxassetnodeindex + 1)));

		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);
			//업데이트후의 값으로 저장하기
			vgrid.back() = (vector<double>(vold, vold + (maxassetnodeindex + 1)));
			ugrid.back() = (vector<double>(uold, uold + (maxassetnodeindex + 1)));
		}
		if (t == vd)
			break;
	}

	double pv_fd;
	pv_fd = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
	//doule check
	auto it_vgrid = vgrid.rbegin();
	auto it_ugrid = ugrid.rbegin();

	pv_fd = intp1d(s0, px, *it_ugrid, 0, maxassetnodeindex, 0);

	cout << "npv_fd " << pv_fd << endl;

	string fn = getFnameTimeStartingWith(string("ts"));
	ofstream fout_ts(fn.c_str());
	fout_ts << "tau,s_tmp,cash,delta,pv,r_for,q,PL" << endl;

	for (long i = 0; i<numMC_; i++)
	{
		vector<double> path;
		vector<double> aPL;
		path.push_back(paras.get_spot());

		auto riter_vgrid = vgrid.rbegin();
		auto riter_ugrid = ugrid.rbegin();

		s_tmp = s0;
		tmpKIFlag = hitFlag;

		int init_i = 0;
		double cash = 0.0;
		double pv = 0.0;
		double delta = 0.0;
		double delta_new = 0.0;
		double PL = 0.0;

		int spot_idx = getIndex(s_tmp, px, 0, maxassetnodeindex);

		if (spot_idx == 0) {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);

			}
			else {
				delta = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}
		}
		else if (spot_idx == maxassetnodeindex) {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
			}
			else {
				delta = ((*riter_ugrid)[spot_idx] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}

		}
		else {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
			}
			else {
				delta = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}
		}

		cash += pv - s_tmp*delta;
		PL = cash - pv + s_tmp*delta;
		aPL.push_back(PL);
		if (0)
			fout_ts << 0 << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				double short_vol = paras.get_Lvol_hybrid(idxT[t - vd], s_tmp);
				double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				for (long t2 = 1; t2 <= daydivide_; t2++) {
					s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
					path.push_back(s_tmp);

					if (s_tmp<kibarrier)
						tmpKIFlag = 1;
				}

				spot_idx = getIndex(s_tmp, px, 0, maxassetnodeindex);

				riter_vgrid++;
				riter_ugrid++;

				if (spot_idx == 0) {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}
				}
				else if (spot_idx == maxassetnodeindex) {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}

				}
				else {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}
				}

				cash *= std::exp(r_forward_p[t - vd] * dt);
				cash -= s_tmp*(delta_new - delta);
				cash += s_tmp*delta*(std::exp(q_forward_p[t - vd] * dt) - 1.0);

				delta = delta_new;
				PL = cash - pv + s_tmp*delta;
				aPL.push_back(PL);
				if (0)
					fout_ts << t - vd << "," << s_tmp << "," << cash << "," << delta << "," << pv << "," << r_forward_p[t - vd] << "," << q_forward_p[t - vd] << "," << PL << endl;
			}

			//autocall payoff 
			double df = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd]);

			if (s_tmp >= autocall_strike[k]) { //check autocallability
											   //mcvalues.push_back(std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]));
				mcvalues.push_back(1.0 + autocall_coupon[k]);

				//pv = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(cash + s_tmp*delta)-mcvalues.back();
				pv = mcvalues.back();

				pvs.push_back(pv);
				PL = cash - pv + s_tmp*delta;
				aPL.back() = PL;
				if (db)
					fout_ts << tmpKIFlag << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;
				break; //k loop
			}

			if (k == nb_autocall) {//we are here because it hasn't been autocalled
				if (s_tmp >= autocall_strike[k]) {
					//mcvalues.push_back(std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]));
					mcvalues.push_back(1.0 + autocall_coupon[k]);

					pv = mcvalues.back();
					pvs.push_back(pv);
					PL = cash - pv + s_tmp*delta;
					aPL.back() = PL;
					if (db)
						fout_ts << tmpKIFlag << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues.push_back(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
						pv = mcvalues.back();
						pvs.push_back(pv);
						PL = cash - pv + s_tmp*delta;
						aPL.back() = PL;
						if (db)
							fout_ts << tmpKIFlag << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

					}
					else if (tmpKIFlag == 0) {
						mcvalues.push_back(1.0 + dummy_coupon);
						pv = mcvalues.back();
						pvs.push_back(pv);
						PL = cash - pv + s_tmp*delta;
						aPL.back() = PL;
						if (db)
							fout_ts << tmpKIFlag << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues.push_back(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					pv = mcvalues.back();
					pvs.push_back(pv);
					PL = cash - pv + s_tmp*delta;
					aPL.back() = PL;
					if (db)
						fout_ts << tmpKIFlag << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;
				}

			} //if k

		}//for k

		paths.push_back(path);
		PLs.push_back(aPL);
		pvs_final_euro.push_back(eop.Simulation3(paras,path,false));
		//eop.Simulation3(paras, path, false);
		if(!db)
			cout << i << endl;
	}//for(i=0..)

	fout_ts.close();

	string fn2 = getFnameTimeStartingWith(string("PL"));
	ofstream fout(fn2.c_str());
	auto it = PLs.begin();
	auto it_euro = pvs_final_euro.begin();
	fout << "numMC,PL" << endl;
	for (auto iter = PLs.begin(); iter != PLs.end(); iter++) {
		fout << iter - it << "," << (*iter).back() << ","<< *it_euro<< endl;
		it_euro++;
	}

	fout.close();

	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] idxT;
	delete[] idxS;
	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
}


void AutocallOption::Simulation3(MarketParameters & paras, std::vector<double>& apath, bool db)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();

	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	vector<vector<double> > paths;

	int hitFlag = this->hitflag;
	double refprice = this->refprice;

	assert(refprice == paras.get_spot());
	assert(paras.get_vdate() == vd);

	assert(paras.get_spot() == apath[0]);

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	vector<double> mcvalues;
	vector<double> pvs;
	vector<vector<double> > PLs;

	int maxassetnodeindex = 400;
	double *px = new double[maxassetnodeindex + 1];
	double *dpx = new double[maxassetnodeindex + 1];
	double *alpha = new double[maxassetnodeindex + 1];
	double *alpha_up = new double[maxassetnodeindex + 1];
	double *alpha_down = new double[maxassetnodeindex + 1];

	double *beta = new double[maxassetnodeindex + 1];
	double *vold = new double[maxassetnodeindex + 1];
	double *vold_up = new double[maxassetnodeindex + 1];
	double *vold_down = new double[maxassetnodeindex + 1];

	double *uold = new double[maxassetnodeindex + 1];
	double *uold_up = new double[maxassetnodeindex + 1];
	double *uold_down = new double[maxassetnodeindex + 1];

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta
	double *uold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *unew = new double[maxassetnodeindex + 1];
	double *unew_up = new double[maxassetnodeindex + 1];
	double *unew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	px[0] = 0.0;
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];

	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	int kiindex = 0;
	kiindex = findlowerindex(px, ThePayoffPtr->GetKiBarrier(), 0, maxassetnodeindex);

	//final BC
	ThePayoffPtr->final_updator(vold, uold, px, 0, maxassetnodeindex);
	for (int i = 0; i <= maxassetnodeindex; i++) {
		vold_up[i] = vold[i];
		vold_down[i] = vold[i];
		uold_up[i] = uold[i];
		uold_down[i] = uold[i];
	}



	/*save vold, uold in vector*/
	vector<vector<double> > vgrid;
	vector<vector<double> > ugrid;

	vgrid.push_back(vector<double>(vold, vold + (maxassetnodeindex + 1)));
	ugrid.push_back(vector<double>(uold, uold + (maxassetnodeindex + 1)));


	if (db) {
		auto it_vgrid = vgrid.rbegin();
		auto it_ugrid = ugrid.rbegin();

		ofstream fuold("uold_final.csv");
		ofstream fvold("vold_final.csv");
		fuold << "PX, uold" << endl;
		fvold << "PX, vold" << endl;

		for (int i = 0; i <= maxassetnodeindex; i++) {
			fuold << px[i] << "," << (*it_ugrid)[i] << endl;
			fvold << px[i] << "," << (*it_vgrid)[i] << endl;
		}

		fuold.close();
		fvold.close();
	}

	signed int t;
	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	int *idxS = new int[maxassetnodeindex + 1];
	for (int i = 0; i <= maxassetnodeindex; i++) {
		idxS[i] = paras.find_index_spot(px[i]);
	}

	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (int k = nb_autocall; k > 0; k--) {
		for (t = autocall_date[k]; t >= std::max(vd + 1, autocall_date[k - 1] + 1); t--) {

			if (t == vd + 1) {
				for (int i = 0; i <= maxassetnodeindex; i++) {
					vold_next[i] = vold[i];  //theta
					uold_next[i] = uold[i];  //theta
				}
			}

			for (int i = 0; i <= maxassetnodeindex; i++) {
				double short_vol = paras.get_Lvol(idxT[t - vd], idxS[i]);
				double short_vol_up = paras.get_Lvol_up(idxT[t - vd], idxS[i]);
				double short_vol_down = paras.get_Lvol_down(idxT[t - vd], idxS[i]);

				alpha[i] = 0.5*short_vol*short_vol*dt;
				alpha_up[i] = 0.5*short_vol_up*short_vol_up*dt;
				alpha_down[i] = 0.5*short_vol_down*short_vol_down*dt;

				beta[i] = (r_forward_p[t - vd] - q_forward_p[t - vd])*dt;
			}

			trimatrix1d(A, B, C, alpha, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_up, B_up, C_up, alpha_up, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);
			trimatrix1d(A_down, B_down, C_down, alpha_down, beta, r_forward_p[t - vd], dt, px, dpx, 1, maxassetnodeindex - 1);

			trimxsolve1d(A, B, C, vold, vnew, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_up, B_up, C_up, vold_up, vnew_up, 0, maxassetnodeindex, 0, 0);
			trimxsolve1d(A_down, B_down, C_down, vold_down, vnew_down, 0, maxassetnodeindex, 0, 0);

			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = vnew[i];
				vold_up[i] = vnew_up[i];
				vold_down[i] = vnew_down[i];
			}

			ThePayoffPtr->copy_v_to_u(vnew, unew, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_up, unew_up, 0, kiindex);
			ThePayoffPtr->copy_v_to_u(vnew_down, unew_down, 0, kiindex);

			trimxsolve1d(A, B, C, uold, unew, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_up, B_up, C_up, uold_up, unew_up, kiindex, maxassetnodeindex, 1, 0);
			trimxsolve1d(A_down, B_down, C_down, uold_down, unew_down, kiindex, maxassetnodeindex, 1, 0);

			ThePayoffPtr->copy_v_to_u(unew, uold, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_down, uold_down, 0, maxassetnodeindex);
			ThePayoffPtr->copy_v_to_u(unew_up, uold_up, 0, maxassetnodeindex);

			vgrid.push_back(vector<double>(vold, vold + (maxassetnodeindex + 1)));
			ugrid.push_back(vector<double>(uold, uold + (maxassetnodeindex + 1)));

		}//for t

		 //update BC
		if (t == autocall_date[k - 1]) {
			if (db) {

				ostringstream oss_uold;
				ostringstream oss_vold;

				oss_uold << "uold" << k-1 << ".csv";
				oss_vold << "vold" << k-1 << ".csv";

				auto it_vgrid = vgrid.rbegin();
				auto it_ugrid = ugrid.rbegin();

				ofstream fuold(oss_uold.str().c_str());
				ofstream fvold(oss_vold.str().c_str());
				fuold << "PX, uold" <<k-1<< endl;
				fvold << "PX, vold" <<k-1<< endl;

				for (int i = 0; i <= maxassetnodeindex; i++) {
					fuold << px[i] << "," << (*it_ugrid)[i] << endl;
					fvold << px[i] << "," << (*it_vgrid)[i] << endl;
				}

				fuold.close();
				fvold.close();
			}
		}

		if (t == autocall_date[k - 1]) {
			ThePayoffPtr->updator(t, vold, uold, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_up, uold_up, px, 0, maxassetnodeindex);
			ThePayoffPtr->updator(t, vold_down, uold_down, px, 0, maxassetnodeindex);
			//업데이트후의 값으로 저장하기
			vgrid.back() = (vector<double>(vold, vold + (maxassetnodeindex + 1)));
			ugrid.back() = (vector<double>(uold, uold + (maxassetnodeindex + 1)));
		}
		if (t == vd)
			break;
	}

	double pv_fd;
	pv_fd = intp1d(s0, px, uold, 1, maxassetnodeindex - 1);
	//doule check
	auto it_vgrid = vgrid.rbegin();
	auto it_ugrid = ugrid.rbegin();

	pv_fd = intp1d(s0, px, *it_ugrid, 0, maxassetnodeindex, 0);

	cout << "npv_fd " << pv_fd << endl;

	ofstream fuold("uold.csv");
	ofstream fvold("vold.csv");
	fuold << "PX, uold" << endl;
	fvold << "PX, vold" << endl;

	for (int i = 0; i <= maxassetnodeindex; i++) {
		fuold << px[i] << "," << (*it_ugrid)[i] << endl;
		fvold << px[i] << "," << (*it_vgrid)[i] << endl;
	}

	fuold.close();
	fvold.close();

	string fn = getFnameTimeStartingWith(string("ts"));
	ofstream fout_ts(fn.c_str());
	fout_ts << "tau,s_tmp,cash,delta,pv,r,q,PL" << endl;

	for (long i = 0; i == 0; i++)
	{
		vector<double> path;
		vector<double> aPL;
		path.push_back(paras.get_spot());

		auto riter_vgrid = vgrid.rbegin();
		auto riter_ugrid = ugrid.rbegin();

		s_tmp = s0;
		tmpKIFlag = hitFlag;

		int init_i = 0;
		double cash = 0.0;
		double pv = 0.0;
		double delta = 0.0;
		double delta_new = 0.0;
		double PL = 0.0;

		int spot_idx = getIndex(s_tmp, px, 0, maxassetnodeindex);

		if (spot_idx == 0) {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);

			}
			else {
				delta = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}
		}
		else if (spot_idx == maxassetnodeindex) {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
			}
			else {
				delta = ((*riter_ugrid)[spot_idx] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}

		}
		else {
			if (tmpKIFlag) {
				delta = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
			}
			else {
				delta = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
				pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
			}
		}

		cash += pv - s_tmp*delta;
		PL = cash - pv + s_tmp*delta;
		aPL.push_back(PL);
		if (db)
			fout_ts << 0 << "," << s_tmp << "," << cash << "," << delta << "," << pv << ",,," << PL << endl;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				s_tmp = apath[t - vd];
				path.push_back(s_tmp);
				if (s_tmp<kibarrier)
					tmpKIFlag = 1;

				spot_idx = getIndex(s_tmp, px, 0, maxassetnodeindex);

				riter_vgrid++;
				riter_ugrid++;

				if (spot_idx == 0) {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx]) / (px[spot_idx + 1] - px[spot_idx]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}
				}
				else if (spot_idx == maxassetnodeindex) {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}

				}
				else {
					if (tmpKIFlag) {
						delta_new = ((*riter_vgrid)[spot_idx + 1] - (*riter_vgrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_vgrid, 0, maxassetnodeindex, 0);
					}
					else {
						delta_new = ((*riter_ugrid)[spot_idx + 1] - (*riter_ugrid)[spot_idx - 1]) / (px[spot_idx + 1] - px[spot_idx - 1]);
						pv = intp1d(s_tmp, px, *riter_ugrid, 0, maxassetnodeindex, 0);
					}
				}

				cash *= std::exp(r_forward_p[t - vd] * dt);
				cash -= s_tmp*(delta_new - delta);
				cash += s_tmp*delta*(std::exp(q_forward_p[t - vd] * dt) - 1.0);

				delta = delta_new;
				PL = cash - pv + s_tmp*delta;
				aPL.push_back(PL);
				if (db)
					fout_ts << t - vd << "," << s_tmp << "," << cash << "," << delta << "," << pv << "," << r_forward_p[t - vd] <<","<< q_forward_p[t - vd] << ","<<  PL << endl;
			}

			//autocall payoff 
			double df = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd]);

			if (s_tmp >= autocall_strike[k]) { //check autocallability
											   //mcvalues.push_back(std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]));
				mcvalues.push_back(1.0 + autocall_coupon[k]);

				//pv = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(cash + s_tmp*delta)-mcvalues.back();
				pv = mcvalues.back();

				pvs.push_back(pv);
				PL = cash - pv + s_tmp*delta;
				aPL.back() = PL;
				if (db)
					fout_ts << i << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << "," << PL << endl;
				break; //k loop
			}

			if (k == nb_autocall) {//we are here because it hasn't been autocalled
				if (s_tmp >= autocall_strike[k]) {
					//mcvalues.push_back(std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]));
					mcvalues.push_back(1.0 + autocall_coupon[k]);

					pv = mcvalues.back();
					pvs.push_back(pv);
					PL = cash - pv + s_tmp*delta;
					aPL.back() = PL;
					if (db)
						fout_ts << i << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << "," << PL << endl;

				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues.push_back(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
						pv = mcvalues.back();
						pvs.push_back(pv);
						PL = cash - pv + s_tmp*delta;
						aPL.back() = PL;
						if (db)
							fout_ts << i << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << "," << PL << endl;

					}
					else if (tmpKIFlag == 0) {
						mcvalues.push_back(1.0 + dummy_coupon);
						pv = mcvalues.back();
						pvs.push_back(pv);
						PL = cash - pv + s_tmp*delta;
						aPL.back() = PL;
						if (db)
							fout_ts << i << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << "," << PL << endl;
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues.push_back(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					pv = mcvalues.back();
					pvs.push_back(pv);
					PL = cash - pv + s_tmp*delta;
					aPL.back() = PL;
					if (db)
						fout_ts << i << "_" << k << "-th Autocalled," << s_tmp << "," << cash << "," << delta << "," << pv << "," << PL << endl;
				}

			} //if k

		}//for k

		paths.push_back(path);
		PLs.push_back(aPL);
	}//for(i=0..)

	fout_ts.close();

	string fn2 = getFnameTimeStartingWith(string("PL"));
	ofstream fout(fn2.c_str());
	auto it = PLs.begin();
	fout << "numMC,PL" << endl;
	for (auto iter = PLs.begin(); iter != PLs.end(); iter++)
		fout << iter - it << "," << (*iter).back() << endl;

	fout.close();

	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] uold;
	delete[] uold_up;
	delete[] uold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;
	delete[] vold_next;
	delete[] unew;
	delete[] unew_up;
	delete[] unew_down;
	delete[] uold_next;

	delete[] A;
	delete[] A_up;
	delete[] A_down;
	delete[] B;
	delete[] B_up;
	delete[] B_down;
	delete[] C;
	delete[] C_up;
	delete[] C_down;

	delete[] idxT;
	delete[] idxS;
	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
}

double AutocallOption::CalcMC_calc2(MarketParam & para, long numMC_)
{
	double s0 = para.get_spot();
	Rate R = para.get_rfrate();
	Rate Q = para.get_q();
	signed int vd = para.get_vdate();
	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	Vol vol = para.get_vol();
	vol.calcLv(s0, R, Q);

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

	//std::random_device rd;
	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	double* mcvalues = new double[numMC_];

	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double* tau_p = new double[autocall_date[nb_autocall]-vd+1];
	double* r_forward_p=new double[autocall_date[nb_autocall]-vd+1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd+1];

	for (signed int i = 0; i <= autocall_date[nb_autocall]-vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = R.getForward(tau_p[i]);
		r_dc_p[i] = R.getIntpRate(tau_p[i]);
		q_forward_p[i] = Q.getForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	for (long i = 0; i<numMC_; i++)
	{
	/*	if (i % (numMC_ / 50) == 0) {
			cout << "now " << i << endl;
		}
*/
		s_tmp = s0;
		tmpKIFlag = hitflag;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				double short_vol = vol.lvol(tau_p[t-vd], s_tmp);

				double drift = (r_forward_p[t-vd] - q_forward_p[t-vd] - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				for (long t2 = 1; t2 <= daydivide_; t2++) {
					s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
					if (s_tmp<kibarrier)
						tmpKIFlag = 1;
				}
			}

			if (s_tmp >= autocall_strike[k]) { //check autocallability
				mcvalues[i] = std::exp(-r_dc_p[autocall_date[k]-vd]*tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]);
				break; //k loop
			}

			//we are here because it is not autocalled at maturity
			if (k == nb_autocall) {

				if (s_tmp >= autocall_strike[k]) {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]);
				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					}
					else if (tmpKIFlag == 0) {
						mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + dummy_coupon);
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
				}

			} //if k

		}//for k

	}//for(i=0..)

	double npv = 0.0;
	for (long i = 0; i<numMC_; i++)
		npv += mcvalues[i];
	npv /= numMC_;

	result[0] = npv;
	result[5] = s0;

	delete[] mcvalues;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;

	return npv;
}

double AutocallOption::CalcMC(MarketParameters & paras, long numMC_)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();
	int nb_autocall = ThePayoffPtr->GetNbAutocall();
	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

	//std::random_device rd;
	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	double* mcvalues = new double[numMC_];

	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];

	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (long i = 0; i<numMC_; i++)
	{
		s_tmp = s0;
		tmpKIFlag = hitflag;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				//double short_vol = paras.lvol(tau_p[t - vd], s_tmp);
				double short_vol = paras.get_Lvol_hybrid(idxT[t - vd], s_tmp);

				double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				for (long t2 = 1; t2 <= daydivide_; t2++) {
					s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
					if (s_tmp<kibarrier)
						tmpKIFlag = 1;
				}
			}

			if (s_tmp >= autocall_strike[k]) { //check autocallability
				mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]);
				break; //k loop
			}

			//we are here because it hasn't been autocalled
			if (k == nb_autocall) {

				if (s_tmp >= autocall_strike[k]) {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]);
				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					}
					else if (tmpKIFlag == 0) {
						mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + dummy_coupon);
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
				}

			} //if k

		}//for k

	}//for(i=0..)

	double npv = 0.0;
	for (long i = 0; i<numMC_; i++)
		npv += mcvalues[i];
	npv /= numMC_;

	result.resize(30, 0.0);
	for (auto iter = result.begin(); iter != result.end(); iter++)
		*iter = 0.0;

	result[0] = npv;
	result[5] = s0;

	delete[] mcvalues;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
	delete[] idxT;
	return npv;
}
double AutocallOption::CalcMC_discrete(MarketParameters & paras, long numMC_)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();
	int nb_autocall = ThePayoffPtr->GetNbAutocall();
	paras.calcLV();

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

	//std::random_device rd;
	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	double* mcvalues = new double[numMC_];

	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double* tau_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* r_dc_p = new double[autocall_date[nb_autocall] - vd + 1];
	double* q_forward_p = new double[autocall_date[nb_autocall] - vd + 1];
	
	double dt = 1 / 365.0;
	for (signed int i = 0; i <= autocall_date[nb_autocall] - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		//q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}
	for (signed int t = autocall_date[nb_autocall]; t >= vd; t--) {
		//CAUTION: t instead of tau
		q_forward_p[t - vd] = paras.getTodayDivAmount(t) / s0 / dt;
	}

	int *idxT = new signed int[autocall_date[nb_autocall] - vd + 1];
	for (int tfv = 0; tfv <= autocall_date[nb_autocall] - vd; tfv++) {
		idxT[tfv] = paras.find_index_term(tfv / 365.0);
	}

	for (long i = 0; i<numMC_; i++)
	{
		s_tmp = s0;
		tmpKIFlag = hitflag;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {

				//double short_vol = paras.lvol(tau_p[t - vd], s_tmp);
				double short_vol = paras.get_Lvol_hybrid(idxT[t - vd], s_tmp);
				double drift = (r_forward_p[t - vd] - q_forward_p[t - vd] - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				for (long t2 = 1; t2 <= daydivide_; t2++) {
					s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
					if (s_tmp<kibarrier)
						tmpKIFlag = 1;
				}
			}

			if (s_tmp >= autocall_strike[k]) { //check autocallability
				mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]);
				break; //k loop
			}

			//we are here because it hasn't been autocalled
			if (k == nb_autocall) {

				if (s_tmp >= autocall_strike[k]) {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + autocall_coupon[k]);
				}
				else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
					}
					else if (tmpKIFlag == 0) {
						mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 + dummy_coupon);
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}
				else {
					mcvalues[i] = std::exp(-r_dc_p[autocall_date[k] - vd] * tau_p[autocall_date[k] - vd])*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
				}

			} //if k

		}//for k

	}//for(i=0..)

	double npv = 0.0;
	for (long i = 0; i<numMC_; i++)
		npv += mcvalues[i];
	npv /= numMC_;

	result.resize(30, 0.0);
	for (auto iter = result.begin(); iter != result.end(); iter++)
		*iter = 0.0;

	result[0] = npv;
	result[5] = s0;

	delete[] mcvalues;

	delete[] tau_p;
	delete[] r_forward_p;
	delete[] r_dc_p;
	delete[] q_forward_p;
	delete[] idxT;
	return npv;
}
double AutocallOption::CalcMC(MarketParam & para, long numMC_)
{
	double s0 = para.get_spot();
	Rate R = para.get_rfrate();
	Rate Q = para.get_q();
	signed int vd = para.get_vdate();
	int nb_autocall = ThePayoffPtr->GetNbAutocall();

	Vol vol = para.get_vol();
	vol.calcLv(s0, R, Q);

	std::vector<signed int> autocall_date;
	autocall_date = ThePayoffPtr->GetAutocall_date();

	double kibarrier = ThePayoffPtr->GetKiBarrier();

	std::vector<double> autocall_strike;
	autocall_strike = ThePayoffPtr->GetAutocall_strike();

	std::vector<double> autocall_coupon;
	autocall_coupon = ThePayoffPtr->GetAutocall_coupon();

//	std::random_device rd;
//	std::mt19937 gen(rd());
	std::mt19937 gen(130);
	std::normal_distribution<>ndist(0, 1);
	//signed int numMC_ = 30000;
	double* mcvalues = new double[numMC_];
	
	double put_strike = ThePayoffPtr->GetPutStrike();
	double dummy_coupon = ThePayoffPtr->GetDummyCoupon();

	double s_tmp;
	int tmpKIFlag;
	int daydivide_ = 1;

	double dt = 1 / 365.0;

	for (long i = 0; i<numMC_; i++)
	{
	/*	if (i % (numMC_ / 50) == 0) {
			cout << "now " << i << endl;
		}
*/
		s_tmp = s0;
	
		tmpKIFlag = hitflag;

		for (int k = 1; k <= nb_autocall; k++) {
			for (signed int t = std::max(autocall_date[k - 1], vd) + 1; t <= autocall_date[k]; t++) {
				double tau = (t - vd) / 365.0;
				

				double short_vol = vol.lvol(tau, s_tmp);
				double r_forward = R.getForward(tau);
				double q_forward = Q.getForward(tau);

				double drift = (r_forward - q_forward - 0.5*short_vol*short_vol)*dt;
				double diff = short_vol*std::sqrt(dt);

				for (long t2 = 1; t2 <= daydivide_; t2++){
					s_tmp = s_tmp*std::exp(drift + diff*ndist(gen));
					if (s_tmp<kibarrier)
						tmpKIFlag = 1;
				}
			}

			if (s_tmp >= autocall_strike[k]){ //check autocallability
				double tau = (autocall_date[k] - vd) / 365.0;
				double dt = 1 / 365.0;

				double r = R.getIntpRate(tau);
				//double r = intp1d((autocalldates_[k] - md) / 365.0, rfrate_term, rfrate, 0, numrfrate - 1);
				//double tau = (autocalldates_[k] - vd) / 365.0;
				mcvalues[i] = std::exp(-r*tau)*(1.0 + autocall_coupon[k]);
				break; //k loop
			}

			//we are here because it is not autocalled at maturity
			if (k == nb_autocall) {
				//on expiry date
				double tau = (autocall_date[k] - vd) / 365.0;
				double dt = 1 / 365.0;
				double r = R.getIntpRate(tau);

				//in case of having knock-in feature
				if (s_tmp >= autocall_strike[k]) {
					mcvalues[i] = std::exp(-r*tau)*(1.0 + autocall_coupon[k]);
				}else if (s_tmp >= kibarrier) {
					if (tmpKIFlag == 1) {
						mcvalues[i] = std::exp(-r*tau)*(1.0 - std::max((put_strike- s_tmp) /refprice, 0.0));
					}
					else if (tmpKIFlag == 0) {
						mcvalues[i] = std::exp(-r*tau)*(1.0 + dummy_coupon);
					}
					else {
						throw std::logic_error("unexpected KIFlag");
					}
				}else {
					mcvalues[i] = std::exp(-r*tau)*(1.0 - std::max((put_strike - s_tmp) / refprice, 0.0));
				}

			} //if k

		}//for k

	}//for(i=0..)

	double npv=0.0;
	for (long i = 0; i<numMC_; i++)
		npv += mcvalues[i];
	npv /= numMC_;

	result[0] = npv; 

	delete[] mcvalues;

	return npv;
}



signed int AutocallOption::GetExpiryd() const
{
	return expiry_date;
}

double AutocallOption::GetRefPrice() const
{
	return refprice;
}

std::vector<double> AutocallOption::GetResult() const
{
	return result;
}

void AutocallOption::PrintResult() const
{
	std::cout.precision(8);
	std::cout << std::fixed;
	cout << "pv " << result[0] << endl;
	cout << "adjdelta " << result[1] << endl;
	cout << "adjgamma " << result[2] << endl;
	cout << "vega " << result[3] << endl;
	cout << "theta " << result[4] << endl;
	cout << "delta(pure) " << result[6] << endl;
	cout << "gamma(pure) " << result[7] << endl;
}

void AutocallOption::PrintSpec() const
{
	cout << "refprice " << refprice << endl;
	cout << "expiry_date " << expiry_date << endl;
	cout << "hitflag " << hitflag << endl;
	ThePayoffPtr->PrintPayoff();
}


double AutocallOption::get_delta(double target, double * px, double* uold, double* vold, int KIFlag, int min, int max, int& init_i)
{
	//search relevant index with init_i
	if (target <= px[min]) {
		init_i = min;
		if (KIFlag == 0) {
			return (uold[min+1]-uold[min]) / (px[min+1]-px[min]);
		}
		else if (KIFlag == 1) {
			return (vold[min + 1] - vold[min]) / (px[min + 1] - px[min]);
		}
	}

	if (target >= px[max]) {
		init_i = max;
		if (KIFlag == 0) {
			return (uold[max] - uold[max-1]) / (px[max] - px[max-1]);
		}
		else if (KIFlag == 1) {
			return (vold[max] - vold[max-1]) / (px[max] - px[max-1]);
		}
	}

	if (px[init_i] <= target && target <px[init_i + 1]) {  //if init_i ==max already, it's returned 
		if (target - px[init_i] <px[init_i + 1] - target) {
			if 
				(KIFlag == 0) {
				return (uold[init_i+1] - uold[init_i-1]) / (px[init_i+1] - px[init_i - 1]);
			}
			else if (KIFlag == 1) {
				return (vold[init_i + 1] - vold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
			}
		}
		else {
			init_i = init_i + 1;
			if (KIFlag == 0) {
				return (uold[init_i + 1] - uold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
			}
			else if (KIFlag == 1) {
				return (vold[init_i + 1] - vold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
			}
		}
	}

	int i = init_i;
	while (1) {   //향후 이부분 개선 
		i = i + 1;
		if (i < max) {
			if (px[i] <= target && target < px[i + 1]) {
				if (target - px[i] < px[i + 1] - target) {
					init_i = i;
				}
				else {
					init_i = i + 1;
				}

				if (KIFlag == 0) {
					return (uold[init_i + 1] - uold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
				}
				else if (KIFlag == 1) {
					return (vold[init_i + 1] - vold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
				}
			}
		}

		int j = init_i - (i - init_i);
		if (j >= 0) {
			if (px[j] <= target && target <px[j + 1]) {
				if (target - px[j] <px[j + 1] - target) {
					init_i = j;
				}
				else {
					init_i = j + 1;
				}

				if (KIFlag == 0) {
					return (uold[init_i + 1] - uold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
				}
				else if (KIFlag == 1) {
					return (vold[init_i + 1] - vold[init_i - 1]) / (px[init_i + 1] - px[init_i - 1]);
				}
			}
		}
	}

	throw std::logic_error("find_index_spot2 - interpolaton fail :findnearestindex, vol strike");
	return -1;
}

double AutocallOption::get_val(double target, double * px, double * uold, double * vold, int KIFlag, int min, int max, int & init_i)
{

	//search relevant index with init_i
	if (target <= px[min]) {
		init_i = min;
		if (KIFlag == 0) {
			return uold[min];
			
		}
		else if (KIFlag == 1) {
			return vold[min];
		}
	}

	if (target >= px[max]) {
		init_i = max;
		if (KIFlag == 0) {
			return uold[max];
		}
		else if (KIFlag == 1) {
			return vold[max];
		}
	}

	if (px[init_i] <= target && target <px[init_i + 1]) {  //if init_i ==max already, it's returned 
		if (KIFlag == 0) {
			return (uold[init_i+1]*(target-px[init_i])+uold[init_i]*(px[init_i+1]-target))/(px[init_i+1]-px[init_i]);
		}
		else if (KIFlag == 1) {
			return (vold[init_i + 1] * (target - px[init_i]) + vold[init_i] * (px[init_i + 1] - target)) / (px[init_i + 1] - px[init_i]);
		}

	}

	int i = init_i;
	while (1) {   //향후 이부분 개선 
		i = i + 1;
		if (i < max) {
			if (px[i] <= target && target < px[i + 1]) {
				init_i = i;

				if (KIFlag == 0) {
					return (uold[init_i + 1] * (target - px[init_i]) + uold[init_i] * (px[init_i + 1] - target)) / (px[init_i + 1] - px[init_i]);
				}
				else if (KIFlag == 1) {
					return (vold[init_i + 1] * (target - px[init_i]) + vold[init_i] * (px[init_i + 1] - target)) / (px[init_i + 1] - px[init_i]);
				}
			}
		}

		int j = init_i - (i - init_i);
		if (j >= 0) {
			if (px[j] <= target && target <px[j + 1]) {
				init_i = j;

				if (KIFlag == 0) {
					return (uold[init_i + 1] * (target - px[init_i]) + uold[init_i] * (px[init_i + 1] - target)) / (px[init_i + 1] - px[init_i]);
				}
				else if (KIFlag == 1) {
					return (vold[init_i + 1] * (target - px[init_i]) + vold[init_i] * (px[init_i + 1] - target)) / (px[init_i + 1] - px[init_i]);
				}
			}
		}
	}

	throw std::logic_error("find_index_spot2 - interpolaton fail :findnearestindex, vol strike");
	return -1;
}

int AutocallOption::getIndex(double target, double * px, int i_min, int i_max) const
{
	if (target <= px[0])
		return (init_i = 0);
	if (target >= px[i_max])
		return (init_i = i_max);

	if (px[init_i] <= target && target <px[init_i + 1]) {
		if (target - px[init_i] <px[init_i + 1] - target) {
			return init_i;
		}
		else {
			return (init_i = init_i + 1);
		}
	}

	int i = init_i;
	while (1) {   //향후 이부분 개선 
		i = i + 1;
		if (i < i_max) {
			if (px[i] <= target && target < px[i + 1]) {
				if (target - px[i] < px[i + 1] - target) {
					return (init_i = i);
				}
				else {
					return (init_i = i + 1);
				}
			}
		}

		int j = init_i - (i - init_i);
		if (j >= 0) {
			if (px[j] <= target && target <px[j + 1]) {
				if (target - px[j] < px[j + 1] - target) {
					return (init_i = j);
				}
				else {
					return (init_i = j + 1);
				}
			}
		}
	}

	throw std::logic_error("find_index_spot2 - interpolaton fail :findnearestindex, vol strike");
	return -1;

}




