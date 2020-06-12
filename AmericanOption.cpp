#include "AmericanOption.h"
#include "j_fd.h"
#include "k_miscellaneous.hpp"
double AmericanOption::Calc(MarketParam& para)
{
	double s0 = para.get_spot();
	
	Rate R = para.get_rfrate();
	Rate Q = para.get_q();

	Vol vol = para.get_vol();
	vol.calcLv(s0, R, Q);

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

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	//int kiindex=0;
	px[0] = 0.0;
	double tmp_ds = refprice*3.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];
	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);

	for (signed int t = expiry_date; t >= para.get_vdate(); t--) {
		if (t == expiry_date) {  //b.c, expiry date
			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = (*ThePayoffPtr)(px[i]);
				vold_up[i] = vold[i];
				vold_down[i] = vold[i];
			}

		}

		if (t == (para.get_vdate() + 1)) {  //for theta, if vd==expiry date ? 
			for (int i = 0; i <= maxassetnodeindex; i++)
				vold_next[i] = vold[i];
		}

		double tau = (t - para.get_vdate()) / 365.0; //time from vdate
		double dt = 1 / 365.0;

		double r_forward = R.getForward(tau);
		double q_forward = Q.getForward(tau);

		for (int i = 0; i <= maxassetnodeindex; i++) {
			double short_vol = vol.lvol(tau, px[i]);
			double short_vol_up = vol.lvol_up(tau, px[i]);
			double short_vol_down = vol.lvol_down(tau, px[i]);

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

		//American feature 
		for (int i = 0; i <= maxassetnodeindex; i++) {
			double po = (*ThePayoffPtr)(px[i]);
			if (vnew[i] <= po) 
				vnew[i] = po;

			if (vnew_up[i] <= po) 
				vnew_up[i] = po;

			if (vnew_down[i] <= po) 
				vnew_down[i] = po;
		}
		//end of american feature

		for (int i = 0; i <= maxassetnodeindex; i++) {
			vold[i] = vnew[i];
			vold_up[i] = vnew_up[i];
			vold_down[i] = vnew_down[i];

		}

	}

	double pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
	double pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);
	double pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
	double pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);

	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;

	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;

	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;

	delete[] vold_next;
	delete[] A;
	delete[] A_up;
	delete[] A_down;

	delete[] B;
	delete[] B_up;
	delete[] B_down;


	delete[] C;
	delete[] C_up;
	delete[] C_down;

	result.resize(30, 0.0);
	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02); //delta
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01); //gamma
	result[3] = 0;  //vega
	result[4] = pv_next - pv;  //theta
	result[5] = s0;
	return pv;

}
double AmericanOption::GetRefPrice() const
{
	return refprice;
}
signed int AmericanOption::GetExpiryd() const
{
	return expiry_date;
}
AmericanOption::AmericanOption(double _refprice, signed int _expiryd, const Payoff& ThePayoff_)
	:refprice(_refprice), expiry_date(_expiryd)
{
	ThePayoffPtr = ThePayoff_.clone();
	result = std::vector<double>(30);
}

std::vector<double> AmericanOption::GetResult() const
{
	return result;
}

AmericanOption::~AmericanOption()
{
	delete ThePayoffPtr;
}

double AmericanOption::Calc(MarketParameters & paras)
{
	//marketparams, conti div
	signed int vd = paras.get_vdate();
	double s0 = paras.get_spot();
	paras.calcLV();

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

	double *vold_next = new double[maxassetnodeindex + 1];  //reserve for theta

	double *vnew = new double[maxassetnodeindex + 1];
	double *vnew_up = new double[maxassetnodeindex + 1];
	double *vnew_down = new double[maxassetnodeindex + 1];

	double *A = new double[maxassetnodeindex + 1];
	double *A_up = new double[maxassetnodeindex + 1];
	double *A_down = new double[maxassetnodeindex + 1];

	double *B = new double[maxassetnodeindex + 1];
	double *B_up = new double[maxassetnodeindex + 1];
	double *B_down = new double[maxassetnodeindex + 1];

	double *C = new double[maxassetnodeindex + 1];
	double *C_up = new double[maxassetnodeindex + 1];
	double *C_down = new double[maxassetnodeindex + 1];

	//int kiindex=0;
	px[0] = 0.0;
	double tmp_ds = refprice*3.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];
	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);
	
	double* tau_p = new double[expiry_date - vd + 1];
	double* r_forward_p = new double[expiry_date - vd + 1];
	double* r_dc_p = new double[expiry_date - vd + 1];
	double* q_forward_p = new double[expiry_date - vd + 1];
	
	for (signed int i = 0; i <= expiry_date - vd; i++) {
		tau_p[i] = (i) / 365.0;
		r_forward_p[i] = paras.getForward(tau_p[i]);
		r_dc_p[i] = paras.getIntpRate(tau_p[i]);
		q_forward_p[i] = paras.getDivForward(tau_p[i]);
	}

	double dt = 1 / 365.0;

	for (signed int t = expiry_date; t >= vd; t--) {
		if (t == expiry_date) {  //b.c, expiry date
			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = (*ThePayoffPtr)(px[i]);
				vold_up[i] = vold[i];
				vold_down[i] = vold[i];
			}
		}

		if (t == (vd + 1)) {  //for theta, if vd==expiry date ? 
			for (int i = 0; i <= maxassetnodeindex; i++)
				vold_next[i] = vold[i];
		}

		for (int i = 0; i <= maxassetnodeindex; i++) {
			double short_vol = paras.lvol(tau_p[t - vd], px[i]);
			double short_vol_up = paras.lvol_up(tau_p[t - vd], px[i]);
			double short_vol_down = paras.lvol_down(tau_p[t - vd], px[i]);

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

		//American feature 
		for (int i = 0; i <= maxassetnodeindex; i++) {
			double po = (*ThePayoffPtr)(px[i]);
			if (vnew[i] <= po)
				vnew[i] = po;

			if (vnew_up[i] <= po)
				vnew_up[i] = po;

			if (vnew_down[i] <= po)
				vnew_down[i] = po;
		}
		//end of american feature

		for (int i = 0; i <= maxassetnodeindex; i++) {
			vold[i] = vnew[i];
			vold_up[i] = vnew_up[i];
			vold_down[i] = vnew_down[i];
		}
	}


	double pv = intp1d(s0, px, vold, 1, maxassetnodeindex - 1);
	double pv_next = intp1d(s0, px, vold_next, 1, maxassetnodeindex - 1);

	//sticky delta, gamma
	double pv_up = intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex - 1);
	double pv_down = intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex - 1);
	//pure delta, gamma
	double pv_up2 = intp1d(s0*1.01, px, vold, 1, maxassetnodeindex - 1);
	double pv_down2 = intp1d(s0*0.99, px, vold, 1, maxassetnodeindex - 1);
	
	delete[] px;
	delete[] dpx;
	delete[] alpha;
	delete[] alpha_up;
	delete[] alpha_down;
	delete[] beta;
	delete[] vold;
	delete[] vold_up;
	delete[] vold_down;
	delete[] vnew;
	delete[] vnew_up;
	delete[] vnew_down;

	delete[] vold_next;
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

	result.resize(30, 0.0);
	result[0] = pv;
	result[1] = (pv_up - pv_down) / (s0*0.02); //sticky delta
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01); //sticky  gamma
	result[3] = 0;  //vega
	result[4] = pv_next - pv;  //theta
	result[5] = 0;
	result[6] = (pv_up2 - pv_down2) / (s0*0.02); //pure delta
	result[7] = (pv_up2 - 2.0*pv + pv_down2) / (s0*0.01) / (s0*0.01); //pure gamma

	return pv;
}
