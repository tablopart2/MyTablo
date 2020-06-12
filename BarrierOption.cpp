#include "BarrierOption.h"
#include "j_fd.h"
#include "k_miscellaneous.hpp"
#include <fstream>
#include <sstream>
#include <ctime>

void print_vold(double *px, double* vold, int maxi)
{
	std::ostringstream oss;

	std::time_t current;
	current = time(NULL);

	oss << std::string("Lvol") << current << std::string(".csv");
	//ossI << std::string("Ivol") << current << std::string(".csv");

	std::ofstream ofsL;
	//std::ofstream ofsI;

	ofsL.open(oss.str().c_str());
	//ofsI.open(ossI.str().c_str());

	for (int i = 0; i<=maxi; i++) {
			ofsL << vold[i] << ",";		
	}

	ofsL.close();
	//ofsI.close();
	oss.clear();
	//ossI.clear();
}

double BarrierOption::Calc(MarketParam& para)
{
	double s0=para.get_spot();
	signed int vd = para.get_vdate();
	
	Rate R=para.get_rfrate();
	Rate Q=para.get_q();
	
	Vol vol=para.get_vol();
	
	vol.calcLv(s0,R,Q);

	int maxassetnodeindex=300;
	double *px=new double[maxassetnodeindex+1];
	double *dpx=new double[maxassetnodeindex+1];
	double *alpha =new double[maxassetnodeindex+1];
	double *alpha_up =new double[maxassetnodeindex+1];
	double *alpha_down =new double[maxassetnodeindex+1];

	double *beta = new double[maxassetnodeindex+1];
	double *vold=new double[maxassetnodeindex+1];
	double *vold_up=new double[maxassetnodeindex+1];
	double *vold_down=new double[maxassetnodeindex+1];

	double *vold_next=new double[maxassetnodeindex+1];  //reserve for theta
	
	double *vnew=new double[maxassetnodeindex+1];
	double *vnew_up=new double[maxassetnodeindex+1];
	double *vnew_down=new double[maxassetnodeindex+1];
	
	double *A=new double[maxassetnodeindex+1];
	double *A_up=new double[maxassetnodeindex+1];
	double *A_down=new double[maxassetnodeindex+1];

	double *B=new double[maxassetnodeindex+1];
	double *B_up=new double[maxassetnodeindex+1];
	double *B_down=new double[maxassetnodeindex+1];

	double *C=new double[maxassetnodeindex+1];
	double *C_up=new double[maxassetnodeindex+1];
	double *C_down=new double[maxassetnodeindex+1];
	
	//int kiindex=0;
	px[0]=0.0;
	double tmp_ds=refprice*2.0/maxassetnodeindex;
	for(int i=1;i<=maxassetnodeindex;i++)
		px[i]=px[i-1]+tmp_ds;
	for(int i=0;i<maxassetnodeindex;i++) //max index of dp is max index of px -1
		dpx[i]=px[i+1]-px[i];
	ThePayoffPtr->ResetFDGrid(px,dpx,1,maxassetnodeindex-1);
	//gridcontrol(px, dpx, 1, maxassetnodeindex-1,strike,0);


	for(signed int t=expiry_date;t>= vd;t--){
		if(t==expiry_date){  //b.c, expiry date
			for(int i=0;i<=maxassetnodeindex;i++){
				//vold[i]=payoff_at_maturity(px[i]);
				vold[i]=(*ThePayoffPtr)(px[i]); 
				vold_up[i]=vold[i];
				vold_down[i]=vold[i];
			}


		}else{
			ThePayoffPtr->updator(vold, px, maxassetnodeindex);
			ThePayoffPtr->updator(vold_up, px, maxassetnodeindex);
			ThePayoffPtr->updator(vold_down, px, maxassetnodeindex);
		}

		if(t==(vd +1)){  //for theta, if vd==expiry date ? 
			for(int i=0;i<=maxassetnodeindex;i++)
					vold_next[i]=vold[i];
		}

		double tau=(t- vd)/365.0; //time from vdate
		double dt=1/365.0;

		double r_forward=R.getForward(tau);
		double q_forward=Q.getForward(tau);

		for(int i=0;i<=maxassetnodeindex;i++){	
			double short_vol=vol.lvol(tau,px[i]);
			double short_vol_up=vol.lvol_up(tau,px[i]);
			double short_vol_down=vol.lvol_down(tau,px[i]);

			alpha[i]=0.5*short_vol*short_vol*dt;
			alpha_up[i]=0.5*short_vol_up*short_vol_up*dt;
			alpha_down[i]=0.5*short_vol_down*short_vol_down*dt;

			beta[i]=(r_forward-q_forward)*dt;
		}

		trimatrix1d(A,B,C,alpha,beta, r_forward, dt, px, dpx,1, maxassetnodeindex-1);
		trimatrix1d(A_up,B_up,C_up,alpha_up,beta, r_forward, dt, px, dpx,1, maxassetnodeindex-1);
		trimatrix1d(A_down,B_down,C_down,alpha_down,beta, r_forward, dt, px, dpx,1, maxassetnodeindex-1);
	
		trimxsolve1d(A,B, C, vold,  vnew, 0, maxassetnodeindex, 0,0);
		trimxsolve1d(A_up,B_up, C_up, vold_up,  vnew_up, 0, maxassetnodeindex, 0,0);
		trimxsolve1d(A_down,B_down, C_down, vold_down,  vnew_down, 0, maxassetnodeindex, 0,0);
		
		for(int i=0;i<=maxassetnodeindex;i++){
			vold[i]=vnew[i];
			vold_up[i]=vnew_up[i];
			vold_down[i]=vnew_down[i];

		}

	}

	
	double pv=intp1d(s0,px,vold,1,maxassetnodeindex-1);
	double pv_next=intp1d(s0,px,vold_next, 1,maxassetnodeindex-1);
	double pv_up=intp1d(s0*1.01, px, vold_up, 1, maxassetnodeindex-1);
	double pv_down=intp1d(s0*0.99, px, vold_down, 1, maxassetnodeindex-1);

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

	result.resize(30,0.0);
	result[0]=pv;
	result[1]=(pv_up-pv_down)/(s0*0.02);
	result[2]=(pv_up-2.0*pv+pv_down)/(s0*0.01)/(s0*0.01);
	result[3]=0;  //vega
	result[4]=pv_next-pv;  //theta
	result[5]=s0;
	return pv;
}

double BarrierOption::Calc(MarketParameters & paras)
{
	double s0 = paras.get_spot();
	signed int vd = paras.get_vdate();
	paras.calcLV();

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
	double tmp_ds = refprice*2.0 / maxassetnodeindex;
	for (int i = 1; i <= maxassetnodeindex; i++)
		px[i] = px[i - 1] + tmp_ds;
	for (int i = 0; i<maxassetnodeindex; i++) //max index of dp is max index of px -1
		dpx[i] = px[i + 1] - px[i];
	ThePayoffPtr->ResetFDGrid(px, dpx, 1, maxassetnodeindex - 1);

	for (signed int t = expiry_date; t >= vd; t--) {
		if (t == expiry_date) {  //b.c, expiry date
			for (int i = 0; i <= maxassetnodeindex; i++) {
				vold[i] = (*ThePayoffPtr)(px[i]);
				vold_up[i] = vold[i];
				vold_down[i] = vold[i];
			}
		}
		else {
			ThePayoffPtr->updator(vold, px, maxassetnodeindex);
			ThePayoffPtr->updator(vold_up, px, maxassetnodeindex);
			ThePayoffPtr->updator(vold_down, px, maxassetnodeindex);
		}

		if (t == (vd + 1)) {  //for theta, if vd==expiry date ? 
			for (int i = 0; i <= maxassetnodeindex; i++)
				vold_next[i] = vold[i];
		}

		double tau = (t - vd) / 365.0; //time from vdate
		double dt = 1 / 365.0;
		double r_forward = paras.getForward(tau);
		double q_forward = paras.getTodayDivAmount(t) / s0 / dt;

		for (int i = 0; i <= maxassetnodeindex; i++) {
			double short_vol = paras.lvol(tau, px[i]);
			double short_vol_up = paras.lvol_up(tau, px[i]);
			double short_vol_down = paras.lvol_down(tau, px[i]);

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
	result[1] = (pv_up - pv_down) / (s0*0.02);
	result[2] = (pv_up - 2.0*pv + pv_down) / (s0*0.01) / (s0*0.01);
	result[3] = 0;  //vega
	result[4] = pv_next - pv;  //theta
	result[5] = s0;

	return pv;
}

double BarrierOption::GetRefPrice() const
{
	return refprice;
}
signed int BarrierOption::GetExpiryd() const
{
	return expiry_date;
}
BarrierOption::BarrierOption(double _refprice, signed int _expiryd,const PayoffBarrier& ThePayoff_)
	:refprice(_refprice),expiry_date(_expiryd)
{
	ThePayoffPtr=ThePayoff_.clone();
	result=std::vector<double>(30);
}

std::vector<double> BarrierOption::GetResult() const
{
	return result;
}

BarrierOption::~BarrierOption()
{
	delete ThePayoffPtr;
}