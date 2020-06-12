
#include <iostream>
#include <fstream> 
#include <string> 
#include <sstream> 
#include <vector>
#include "volatility.h"
#include "EtcFunctions.h"
#include <cmath> 
#include <ctime> 


void Vol::foutLvol() const
{
	std::ostringstream ossL;
	std::ostringstream ossI;
	
	std::time_t current;
	current=time(NULL);

	ossL<<std::string("Lvol")<<current<<std::string(".csv");
	ossI<<std::string("Ivol")<<current<<std::string(".csv");

	std::ofstream ofsL;
	std::ofstream ofsI;
	
	ofsL.open(ossL.str().c_str());
	ofsI.open(ossI.str().c_str());
	for(int i=0;i<nb_vol_term;i++){
		ofsL<<"\n";
		ofsI<<"\n";
		for(int j=0;j<nb_vol_strike;j++){
			ofsL<<Lvol[i][j]<<",";
			ofsI<<Ivol[i][j]<<",";

		}
	}
	ofsL.close();
	ofsI.close();
	ossL.clear();
	ossI.clear();

}
void Vol::set_const_vol(double v)
{
	for(int i=0;i<nb_vol_term;i++)
		for(int j=0;j<nb_vol_strike;j++)
			Ivol[i][j]=v;
}

int Vol::find_index_term(double target) const
{

	if(target<=vol_term[0])
		return 0;
	if(target>=vol_term[nb_vol_term-1])
		return nb_vol_term-1;

	for(int i=0;i<nb_vol_term-1;i++){
		if(vol_term[i]<=target && target <vol_term[i+1]){
			if(target-vol_term[i] <vol_term[i+1]-target){
				return i;
			}else{
				return i+1;
			}
		}
	}

	throw std::logic_error("interpolaton fail :findnearestindex, vol term");
	return -1;
}

int Vol::find_index_spot(double target) const
{
	
	if(target<=vol_strike[0])
		return 0;
	if(target>=vol_strike[nb_vol_strike-1])
		return nb_vol_strike-1;

	//if (vol_strike[i_s] <= target && target <vol_strike[i_s + 1]) {
	//	if (target - vol_strike[i_s] <vol_strike[i_s + 1] - target) {
	//		return i_s;
	//	}
	//	else {
	//		return i_s + 1;
	//	}
	//}

	for(int i=0;i<nb_vol_strike-1;i++){
		if(vol_strike[i]<=target && target <vol_strike[i+1]){
			if(target-vol_strike[i] <vol_strike[i+1]-target){
				return i;
			}else{
				return i+1;
			}
		}
	}

	throw std::logic_error("interpolaton fail :findnearestindex, vol strike");
	return -1;
}

int Vol::find_index_spot2(double target) const
{
	if (target <= vol_strike[0])
		return (init_i=0);
	if (target >= vol_strike[nb_vol_strike - 1])
		return (init_i=nb_vol_strike - 1);

	if (vol_strike[init_i] <= target && target <vol_strike[init_i + 1]) {  // out of range??
		if (target - vol_strike[init_i] <vol_strike[init_i + 1] - target) {
			return init_i;
		}
		else {
			return (init_i = init_i + 1);
		}
	}

	int i = init_i;
	while (1) {   //향후 이부분 개선 
		i = i + 1;
		if (i < nb_vol_strike-1) {
			if (vol_strike[i] <= target && target < vol_strike[i + 1]) {
				if (target - vol_strike[i] < vol_strike[i + 1] - target) {
					return (init_i = i);
				}
				else {
					return (init_i = i + 1);
				}
			}
		}

		int j = init_i - (i - init_i);
		if (j >= 0) {
			if (vol_strike[j] <= target && target <vol_strike[j + 1]) {
				if (target - vol_strike[j] <vol_strike[j + 1] - target) {
					return (init_i = j);
				}
				else {
					return (init_i = j + 1);
				}
			}
		}
	}

	//for (int i = init_i + 1; i<nb_vol_strike - 1; i++) {
	//	if (vol_strike[i] <= target && target <vol_strike[i + 1]) {
	//		if (target - vol_strike[i] <vol_strike[i + 1] - target) {
	//			return (init_i = i);
	//		}
	//		else {
	//			return (init_i = i + 1);
	//		}
	//	}

	//	int j = init_i - (i - init_i);
	//	if (j > 0) {
	//		if (vol_strike[j] <= target && target <vol_strike[j + 1]) {
	//			if (target - vol_strike[j] <vol_strike[j + 1] - target) {
	//				return (init_i = j);
	//			}
	//			else {
	//				return (init_i = j + 1);
	//			}
	//		}
	//	}
	//}



	throw std::logic_error("find_index_spot2 - interpolaton fail :findnearestindex, vol strike");
	return -1;
}




void Vol::calcLv(double spot, const Rate& R, const Rate& Q)
{
	EtcFunctions fct;
	std::vector<std::vector<double> >dvol_dt(nb_vol_term,std::vector<double>(nb_vol_strike));
	std::vector<std::vector<double> >dvol_dstrike(nb_vol_term,std::vector<double>(nb_vol_strike));
	std::vector<std::vector<double> >dvol2_dstrike2(nb_vol_term,std::vector<double>(nb_vol_strike));

	double spot_up=spot*1.01;
	double spot_down=spot*0.99;

    for(int i=0;i<nb_vol_term;i++){
        if(i==0){
            for(int j=0;j<nb_vol_strike;j++)
                dvol_dt[i][j]=(Ivol[i+1][j]-Ivol[i][j])/(vol_term[i+1]-vol_term[i]);
        }else if(i==(nb_vol_term-1)){
            for(int j=0;j<nb_vol_strike;j++)
                dvol_dt[i][j]=(Ivol[i][j]-Ivol[i-1][j])/(vol_term[i]-vol_term[i-1]);
        }else{
            for(int j=0;j<nb_vol_strike;j++)
                dvol_dt[i][j]=(Ivol[i+1][j]-Ivol[i-1][j])/(vol_term[i+1]-vol_term[i-1]);
        }
    }

   for(int j=0;j<nb_vol_strike;j++){
        if(j==0){
            for(int i=0;i<nb_vol_term;i++){
                dvol_dstrike[i][j]=(Ivol[i][j+1]-Ivol[i][j])/(vol_strike[j+1]-vol_strike[j]);
                dvol2_dstrike2[i][j]=(Ivol[i][j+2]-2.0*Ivol[i][j+1]+Ivol[i][j])/(vol_strike[j+2]-vol_strike[j+1])/(vol_strike[j+1]-vol_strike[j]);
            }
        }else if(j==nb_vol_strike-1){
            for(int i=0;i<nb_vol_term;i++){
                dvol_dstrike[i][j]=(Ivol[i][j]-Ivol[i][j-1])/(vol_strike[j]-vol_strike[j-1]);
                dvol2_dstrike2[i][j]=(Ivol[i][j]-2.0*Ivol[i][j-1]+Ivol[i][j-2])/(vol_strike[j]-vol_strike[j-1])/(vol_strike[j-1]-vol_strike[j-2]);
            }
        }else{
            for(int i=0;i<nb_vol_term;i++){
                dvol_dstrike[i][j]=(Ivol[i][j+1]-Ivol[i][j-1])/(vol_strike[j+1]-vol_strike[j-1]);
                dvol2_dstrike2[i][j]=(Ivol[i][j+1]-2.0*Ivol[i][j]+Ivol[i][j-1])/(vol_strike[j+1]-vol_strike[j])/(vol_strike[j]-vol_strike[j-1]);
            }
        }
    }

    for(int i=0;i<nb_vol_term;i++)
        for(int j=0;j<nb_vol_strike;j++){
            double v, dvdk,dvdt,dv2dk2;
            double d1,d2,d1_up,d1_down, d2_up,  d2_down,t,k,r,q;
            v=Ivol[i][j];
            dvdk=dvol_dstrike[i][j];
            dvdt=dvol_dt[i][j];
            dv2dk2=dvol2_dstrike2[i][j];
            t=vol_term[i];
            k=vol_strike[j];
			r=R.getIntpRate(t);
			q=Q.getIntpRate(t);
			//std::cout << std::sqrt(5.0);
			
            //r=interpolation1D(r_term,r_curve,nb_r,t);
            //q=interpolation1D(div_term,div_curve,nb_div,t);
			
            d1=(log(spot/k)+(r-q+0.5*v*v)*t)/(v*std::sqrt(t));
			d1_up=(log(spot_up/k)+(r-q+0.5*v*v)*t)/(v*std::sqrt(t));
			d1_down=(log(spot_down/k)+(r-q+0.5*v*v)*t)/(v*std::sqrt(t));

            d2=d1-v*std::sqrt(t);
            d2_up=d1_up-v*std::sqrt(t);
			d2_down=d1_down-v*std::sqrt(t);
			Lvol[i][j]=std::sqrt(v*v+2.0*t*v*dvdt+2*(r-q)*k*t*v*dvdk)/std::sqrt(1.0+2.0*k*d1*std::sqrt(t)*dvdk+k*k*t*(d1*d2*dvdk*dvdk+v*dv2dk2));
			Lvol_up[i][j]=std::sqrt(v*v+2.0*t*v*dvdt+2*(r-q)*k*t*v*dvdk)/std::sqrt(1.0+2.0*k*d1_up*std::sqrt(t)*dvdk+k*k*t*(d1_up*d2_up*dvdk*dvdk+v*dv2dk2));
			Lvol_down[i][j]=std::sqrt(v*v+2.0*t*v*dvdt+2*(r-q)*k*t*v*dvdk)/std::sqrt(1.0+2.0*k*d1_down*std::sqrt(t)*dvdk+k*k*t*(d1_down*d2_down*dvdk*dvdk+v*dv2dk2));
        }

}

double Vol::lvol_down(double tau, double s) const
{
//interpolate local vol from Lvol surface
	int i_t=find_index_term(tau);
	int i_s=find_index_spot(s);

	return Lvol_down[i_t][i_s];
}
double Vol::lvol_up(double tau, double s) const
{
//interpolate local vol from Lvol surface
	int i_t=find_index_term(tau);
	int i_s=find_index_spot(s);

	return Lvol_up[i_t][i_s];
}
double Vol::lvol(double tau, double s) const
{
//(flat) interpolate local vol from Lvol surface
	int i_t=find_index_term(tau);
	int i_s=find_index_spot(s);

	return Lvol[i_t][i_s];
}

double Vol::lvol2(double t, double s)
{
	return 0.0;
}

double Vol::getInpVol(double t, double k) const
{//(flat) interpolation from Ivol surface
	int i_t=find_index_term(t);
	int i_k=find_index_spot(k);
	return Ivol[i_t][i_k];
}

double Vol::getBSVol(double t, double k) const
{//(linear) interpolation from Ivol surface
	if(t<=vol_term[0]){
		t=vol_term[0];
	}else if(t>=vol_term[nb_vol_term-1]){
		t=vol_term[nb_vol_term-1];
	}

	if(k<=vol_strike[0]){
		k=vol_strike[0];
	}else if(k>=vol_strike[nb_vol_strike-1]){
		k=vol_strike[nb_vol_strike-1];
	}

	for(int i=1;i<nb_vol_term;i++)
		for(int j=1;j<nb_vol_strike;j++)
			if(t<=vol_term[i] && k<=vol_strike[j]){
				double v=Ivol[i][j]*(t-vol_term[i-1])*(k-vol_strike[j-1])+
					Ivol[i][j-1]*(t-vol_term[i-1])*(vol_strike[j]-k)+
					Ivol[i-1][j]*(vol_term[i]-t)*(k-vol_strike[j-1])+
					Ivol[i-1][j-1]*(vol_term[i]-t)*(vol_strike[j]-k);
				return v/(vol_strike[j]-vol_strike[j-1])/(vol_term[i]-vol_term[i-1]);
			}

	throw std::logic_error("linear vol interpolaton failed, may come from bugs");
}

void Vol::set_vol_term(double* v, int n)
{
	if(n!=nb_vol_term)
		throw;
	for(int i=0;i<n;i++)
		vol_term[i]=v[i];
		//vol_term.push_back(v[i]);

}

void Vol::set_vol_strike(double* v, int n)
{
	if(n!=nb_vol_strike)
		throw;
	for(int i=0;i<n;i++)
		//vol_strike.push_back(v[i]);
		vol_strike[i]=v[i];
}

void Vol::print() const
{
	using namespace std;
	cout<<fixed;
	cout.precision(4);
	
	for(int i=0;i<nb_vol_term;i++){
		cout<<endl;
		for(int j=0;j<nb_vol_strike;j++)
			cout<<Ivol[i][j]<<" ";
	}
}

void Vol::printLV() const
{
	using namespace std;
	cout<<fixed;
	cout.precision(4);
	
	for(int i=0;i<nb_vol_term;i++){
		cout<<endl;
		for(int j=0;j<nb_vol_strike;j++)
			cout<<Lvol[i][j]<<" ";
	}
}

void Vol::set_volsurface_by_term(double* row,int row_num, int nb_strike)
{
//	Ivol.push_back(std::vector<double>());
	for(int i=0;i<nb_strike;i++)
		//Ivol.back().push_back(row[i]);
		Ivol[row_num][i]=row[i];

}

void Vol::set_vol_by_point(double* arr, int nb_term, int nb_strike)
{
	for(int i=0;i<nb_term;i++)
		for(int j=0;j<nb_strike;j++)
			Ivol[i][j]=*(arr+i+j*nb_term);  //CHECK AGAIN, COMPARE  WITH VBA VERSION
}
void Vol::set_vol_by_point_vba(double * arr, int nb_term, int nb_strike)
{
	//this initializer is for VB-sytle(column-based) arrray 
	for (int i = 0; i<nb_term; i++)
		for (int j = 0; j<nb_strike; j++)
			Ivol[i][j] = *(arr + i + j*nb_term);
}
void Vol::set_vol_by_point_vba2(double * arr, int nb_term, int nb_strike)
{
	//this initializer is for VB-sytle(column-based) arrray 
	//In case arr starts from index 0 ,modified 2019.01.23 
	//*(arr+12)==[0][0], *(arr+21)==[9][0], *(arr+32)==[9][1]
	for (int i = 0; i<nb_term; i++)
		for (int j = 0; j<nb_strike; j++)
			Ivol[i][j] = *(arr + (i + 1) + (j + 1)*(nb_term+1));
}

Vol::Vol(int _nb_vol_term, int _nb_vol_strike):nb_vol_term(_nb_vol_term),nb_vol_strike(_nb_vol_strike)
{
	//Ivol=std::vector<<std::vector<double> >(nb_vol_term, std::vector<double>(nb_vol_strike));
	/*for(int i=0;i<nb_vol_term;i++)
		Ivol.push_back(std::vector<double>(nb_vol_strike,)
*/
	//std::vector<std::vector<double> > Ivol(nb_vol_strike,std::vector<double>(nb_vol_term,0.0));
	//std::vector<double> vol_strike(nb_vol_strike);
	//std::vector<double> vol_term(nb_vol_term);
	vol_strike=std::vector<double>(nb_vol_strike);
	vol_term=std::vector<double>(nb_vol_term);
	Ivol=std::vector<std::vector<double> >(nb_vol_term, std::vector<double>(nb_vol_strike));
	Lvol=std::vector<std::vector<double> >(nb_vol_term, std::vector<double>(nb_vol_strike));
	Lvol_up=std::vector<std::vector<double> >(nb_vol_term, std::vector<double>(nb_vol_strike));
	Lvol_down=std::vector<std::vector<double> >(nb_vol_term, std::vector<double>(nb_vol_strike));
	
	//	std::vector<std::vector<double> > Lvol(nb_vol_strike,std::vector<double>(nb_vol_term,0.0));

	//i_t = 0;
	//i_s = 0;
}

Vol::Vol(vector<double>& vol_term_, vector<double>& vol_strike_, vector<vector<double> >& Ivol_):vol_strike(vol_strike_),vol_term(vol_term_), Ivol(Ivol_)
{
	Lvol = Ivol;
	Lvol_up = Ivol;
	Lvol_down = Ivol;
	nb_vol_term = vol_term_.size();;
	nb_vol_strike = vol_strike_.size();
}

Vol::~Vol()
{
	
}

void Vol::Vol_up(double bumping)
{
	//Vol vol(nb_vol_term, nb_vol_strike);

	for(int i=0;i<nb_vol_term;i++)
		for(int j=0;j<nb_vol_strike;j++)
			Ivol[i][j]+=0.01;
}

void Vol::set_Ivol(int i_t, int i_k, double v)
{
	if(i_t>=nb_vol_term || i_k >=nb_vol_strike) 
		throw std::logic_error("out of range, Vol::setIvol");
	Ivol[i_t][i_k]=v;
}

void Vol::reset_Ivol_up()
{
	for (int i = 0; i < nb_vol_term;i++)
		for (int j = 0; j < nb_vol_strike;j++)
			Ivol[i][j] += 0.01;
}

double Vol::get_Ivol(int i_t, int i_k) const
{
	if(i_t>=nb_vol_term || i_k >=nb_vol_strike) 
		throw std::logic_error("out of range, Vol::setIvol");
	return Ivol[i_t][i_k];
}

double Vol::get_Lvol(int i_t, int i_k) const
{
	if (i_t >= nb_vol_term || i_k >= nb_vol_strike)
		throw std::logic_error("out of range, Vol::setIvol");
	return Lvol[i_t][i_k];
}

double Vol::get_Lvol_up(int i_t, int i_k) const
{
	if (i_t >= nb_vol_term || i_k >= nb_vol_strike)
		throw std::logic_error("out of range, Vol::setIvol");
	return Lvol_up[i_t][i_k];
}

double Vol::get_Lvol_down(int i_t, int i_k) const
{
	if (i_t >= nb_vol_term || i_k >= nb_vol_strike)
		throw std::logic_error("out of range, Vol::setIvol");
	return Lvol_down[i_t][i_k];
}

//void Vol::setVol(char* fname)
//{
//	std::ifstream fs(fname);  //filestream
//	if(fs.bad())
//		throw std::domain_error("no file");
//	std::string row;//line segment
//	std::string v;//voldata
//	long i,j;
//
//	std::getline(fs,row);//first row is k/s
//	std::istringstream iss(row);
//	std::getline(iss,v,','); //eliminate upper left charactor, here, t\k
//
//	j=0;
//	while(std::getline(iss,v,',')){
//		//volmoneyness[j]=atof(v.c_str());
//		j++;
//	}
//
//	i=0; //line number
//	while(std::getline(fs,row)){
//		std::istringstream iss(row);
//		std::getline(iss,v,','); //v is the first element of each row, tenor
//		//volterm[i]=atof(v.c_str());
//		j=0;//data number in the line
//		while(std::getline(iss,v,','))
//		{
//			//voldata[nummoneyness*i+j]=atof(v.c_str())/100.0;
//			j++;
//		}
//		i++;
//	}
//}
//
//void ivlv(double spot, double **Lvol_kospi,
//               double **Ivol_kospi,
//               double* vol_term, int nb_vol_term, double* vol_strike, int nb_vol_strike,
//               double* r_curve, double* r_term,int nb_r,
//               double* div_curve, double* div_term, int nb_div)
//{
////    myfunction fct;
//    double** dvol_dt;
//    double** dvol_dstrike;
//    double** dvol2_dstrike2;
//    dvol_dt=new double*[nb_vol_term];
//    dvol_dstrike=new double*[nb_vol_term];
//    dvol2_dstrike2=new double*[nb_vol_term];
//
//    for(int i=0;i<nb_vol_term;i++){
//        dvol_dt[i]=new double[nb_vol_strike];
//        dvol_dstrike[i]=new double[nb_vol_strike];
//        dvol2_dstrike2[i]=new double[nb_vol_strike];
//    }
//
//
////    double dvol_dt[10][17];
////    double dvol_dstrike[10][17];
////    double dvol2_dstrike2[10][17];
//    //dv/dt
//    for(int i=0;i<nb_vol_term;i++){
//        if(i==0){
//            for(int j=0;j<nb_vol_strike;j++)
//                dvol_dt[i][j]=(Ivol_kospi[i+1][j]-Ivol_kospi[i][j])/(vol_term[i+1]-vol_term[i]);
//        }else if(i==(nb_vol_term-1)){
//            for(int j=0;j<nb_vol_strike;j++)
//                dvol_dt[i][j]=(Ivol_kospi[i][j]-Ivol_kospi[i-1][j])/(vol_term[i]-vol_term[i-1]);
//        }else{
//            for(int j=0;j<nb_vol_strike;j++)
//                dvol_dt[i][j]=(Ivol_kospi[i+1][j]-Ivol_kospi[i-1][j])/(vol_term[i+1]-vol_term[i-1]);
//        }
//    }
//
//    //dv/dk
//    for(int j=0;j<nb_vol_strike;j++){
//        if(j==0){
//            for(int i=0;i<nb_vol_term;i++){
//                dvol_dstrike[i][j]=(Ivol_kospi[i][j+1]-Ivol_kospi[i][j])/(vol_strike[j+1]-vol_strike[j]);
//                dvol2_dstrike2[i][j]=(Ivol_kospi[i][j+2]-2.0*Ivol_kospi[i][j+1]+Ivol_kospi[i][j])/(vol_strike[j+2]-vol_strike[j+1])/(vol_strike[j+1]-vol_strike[j]);
//            }
//        }else if(j==nb_vol_strike-1){
//            for(int i=0;i<nb_vol_term;i++){
//                dvol_dstrike[i][j]=(Ivol_kospi[i][j]-Ivol_kospi[i][j-1])/(vol_strike[j]-vol_strike[j-1]);
//                dvol2_dstrike2[i][j]=(Ivol_kospi[i][j]-2.0*Ivol_kospi[i][j-1]+Ivol_kospi[i][j-2])/(vol_strike[j]-vol_strike[j-1])/(vol_strike[j-1]-vol_strike[j-2]);
//            }
//        }else{
//            for(int i=0;i<nb_vol_term;i++){
//                dvol_dstrike[i][j]=(Ivol_kospi[i][j+1]-Ivol_kospi[i][j-1])/(vol_strike[j+1]-vol_strike[j-1]);
//                dvol2_dstrike2[i][j]=(Ivol_kospi[i][j+1]-2.0*Ivol_kospi[i][j]+Ivol_kospi[i][j-1])/(vol_strike[j+1]-vol_strike[j])/(vol_strike[j]-vol_strike[j-1]);
//            }
//        }
//    }
//
//    for(int i=0;i<nb_vol_term;i++)
//        for(int j=0;j<nb_vol_strike;j++){
//            double v, dvdk,dvdt,dv2dk2;
//            double d1,d2,t,k,r,q;
//            v=Ivol_kospi[i][j];
//            dvdk=dvol_dstrike[i][j];
//            dvdt=dvol_dt[i][j];
//            dv2dk2=dvol2_dstrike2[i][j];
//            t=vol_term[i];
//            k=vol_strike[j];
//            //r=interpolation1D(r_term,r_curve,nb_r,t);
//            //q=interpolation1D(div_term,div_curve,nb_div,t);
//            d1=(log(spot/k)+(r-q+0.5*v*v)*t)/(v*sqrt(t));
//            d2=d1-v*sqrt(t);
//            Lvol_kospi[i][j]=sqrt(v*v+2.0*t*v*dvdt+2*(r-q)*k*t*v*dvdk)/sqrt(1.0+2.0*k*d1*sqrt(t)*dvdk+k*k*t*(d1*d2*dvdk*dvdk+v*dv2dk2));
//
//        }
//
//    for(int i=0;i<nb_vol_term;i++){
//        delete[] dvol_dt[i];
//        delete[] dvol_dstrike[i];
//        delete[] dvol2_dstrike2[i];
//    }
//
//    delete[] dvol_dt;
//    delete[] dvol_dstrike;
//    delete[] dvol2_dstrike2;
//
//
//
//}