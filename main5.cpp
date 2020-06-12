#include <iostream>
#include <cstdlib>
#include <cassert>
#include <random>
#include <vector> 
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>
#include <string>

#include "volatility.h"
#include "MarketParam.h"
#include "MarketParameters.h"
#include "PayoffPut.h"
#include "PayoffCall.h"
#include "EuropeanOption.h"
#include "AmericanOption.h"
#include "k_OptionFormula.hpp"
#include "k_miscellaneous.hpp"
#include "PayoffBarCUO.h"
#include "BarrierOption.h"
#include "AutocallOption.h"
#include "AsianOption.h"
#include "Dividend.h"
#include "j_fd.h"

using namespace std;

enum FvsS{flat,surface};
enum CalcMode{calc1,calc2};

MarketParameters init_paras_file(const char* vol_csv, const char* rate_csv, const char* div_csv);
MarketParameters init_paras_file(const char* vol_csv, const char* rate_csv, const char* div_csv, const char* div_discrete_csv);

void test_american(double spot)
{

	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.54,1.64,1.72,1.77,1.8,1.85,1.88,1.90,1.92,1.92,1.92,1.92,1.92,1.92,1.92 };
	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.0191780821917808,0.038356164,0.057534247,0.084931507,0.104109589,0.123287671,0.142465753,0.167123288,0.252054795,0.334246575,0.419178082,0.504109589,0.580821918,0.665753425,0.747945205,0.832876712,0.915068493,1,1.084931507,1.167123288,1.252054795,1.334246575,1.419178082,1.504109589,1.583561644,1.668493151,1.750684932,1.835616438,1.917808219,	2.002739726	,2.087671233,2.169863014,2.254794521,2.336986301,2.421917808,2.506849315,2.583561644,2.668493151,2.750684932,2.835616438,2.917808219,3.002739726 };
	double arr_q[] = { 0.00000,0.00000,0.00000,2.57097, 2.09737,1.77112,1.53270, 1.30656 ,0.86631 ,5.16283, 4.11677,3.42318,3.34816,2.92103,2.60139,2.67302,2.43293,2.22630,2.26258, 2.10325,1.96057,3.00733,2.82806 ,2.66837, 2.53449,2.53675,2.41823 ,2.46233 ,2.35680 ,2.25685 ,2.27446 ,2.18831 ,2.10588 ,2.66917,2.57598,2.48871,2.41481,2.42319,2.35115 ,2.38170 ,2.31461 ,2.24914 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08493150,0.16712330,0.25205480,0.50410960,0.74794520,1.00000000,1.50410960,2.00273970,2.50684930,3.00273970 }; //10terms
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3095,0.2907,0.2719,0.2528,0.2331,0.2136,0.1816,0.1388,0.1081,0.105,0.1296,0.1606,0.1795,0.1892,0.1978,0.2057,0.2129 };
	double t1[] = { 0.2807,0.2619,0.2431,0.2239,0.2042,0.1845,0.1586,0.1285,0.1095,0.1027,0.1076,0.1306,0.1518,0.1634,0.172,0.1798,0.1869 };
	double t2[] = { 0.2621,0.2448,0.2274,0.2098,0.1916,0.172,0.1505,0.128,0.1131,0.1076,0.1125,0.1264,0.1413,0.1517,0.16,0.1676,0.1745 };
	double t3[] = { 0.2376,0.2219,0.2061,0.1902,0.1741,0.1576,0.141,0.1261,0.1169,0.1161,0.1215,0.1292,0.1371,0.1447,0.1517,0.1583,0.1644 };
	double t4[] = { 0.2286,0.214,0.1997,0.1854,0.1712,0.1571,0.1437,0.1319,0.1238,0.1211,0.1231,0.128,0.1339,0.14,0.146,0.1517,0.1571 };
	double t5[] = { 0.221,0.2077,0.1946,0.1819,0.1694,0.1574,0.1462,0.1366,0.1295,0.1261,0.1261,0.1287,0.1327,0.1374,0.1422,0.1469,0.1516 };
	double t6[] = { 0.2063,0.1949,0.1841,0.1737,0.1639,0.1549,0.1469,0.1403,0.1355,0.1329,0.1322,0.1331,0.1352,0.1379,0.1409,0.1442,0.1475 };
	double t7[] = { 0.1993,0.1897,0.1806,0.1721,0.1642,0.157,0.1506,0.1453,0.1411,0.1383,0.1368,0.1364,0.1369,0.1381,0.1398,0.1418,0.144 };
	double t8[] = { 0.1927,0.1844,0.1768,0.1696,0.1631,0.1572,0.1521,0.1477,0.1443,0.1418,0.1402,0.1394,0.1393,0.1398,0.1407,0.1419,0.1434 };
	double t9[] = { 0.1892,0.182,0.1754,0.1692,0.1637,0.1586,0.1542,0.1505,0.1474,0.145,0.1432,0.1421,0.1415,0.1414,0.1417,0.1423,0.1431 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);

	double refprice = 297.22;
	double rf_for_bs = r.getIntpRate(1.0);
	double div_for_bs = q.getIntpRate(1.0);
	double putstrike = 297.22*0.9;
	double callstrike = 297.22*1.1;
	double intpvol = volat.getInpVol(1.0, putstrike);

	MarketParam para(0, spot, volat, r, q);

	PayoffPut putpay(putstrike);
	//EuropeanOption EurPut(refprice, 365, putpay);
	//EurPut.Calc(para);
	AmericanOption AmPut(refprice, 365, putpay);
	AmPut.Calc(para);

	/*EuropeanVanilla Vput(putpay,refprice, ed);

	EuropeanPut eput(refprice,putstrike,365);
	EuropeanCall ecall(refprice,callstrike,365);

	*/
	std::cout.precision(8);
	std::cout << std::fixed;

	//std::cout <<"eput:  s=" <<spot<<" "<< eput.pv(para)<<" vol " << volat.getInpVol(365/365,putstrike)<< "r " << r.getIntpRate(1.0) << " q " << q.getIntpRate(1.0) <<  std::endl;
	//std::cout <<"ecall: "<< ecall.pv(para)<< "s= " <<spot<<" "<<" vol " <<volat.getInpVol(365/365,callstrike)<< " r " << r.getIntpRate(1.0) << " q " << q.getIntpRate(1.0) << std::endl;
	std::vector<double> rs = AmPut.GetResult();

	std::cout << "pv(AMERICAN) " << rs[0] << std::endl;
	std::cout << "delta " << rs[1] << std::endl;
	std::cout << "gamma " << rs[2] << std::endl;
	std::cout << "spot " << rs[5] << std::endl;
	std::cout << "BS vol=" << volat.getBSVol(365 / 365, putpay.GetStrike()) << std::endl;
	std::cout << "implied vol=" << put_iv(spot, putstrike, rf_for_bs, div_for_bs, 0.1, 1.0, rs[0]) << std::endl;
	//std::cout << "delta(closed, Flatvol)=" << k_fmla_BSPutDelta(spot, putstrike, rf_for_bs, div_for_bs, volat.getBSVol(1.0, putpay.GetStrike()), 1.0) << std::endl;
	//std::cout << "gamma(closed, Flatvol)=" << k_fmla_BSGamma(spot, putstrike, rf_for_bs, div_for_bs, volat.getBSVol(1.0, putpay.GetStrike()), 1.0) << std::endl;

	std::cout << "\n";

	//Vol volat_up(volat);
	//volat_up.Vol_up(0.01);
	//MarketParam para_volup(0, spot, volat_up, r, q);
	////EuropeanOption EurPut_vega(refprice, 365, putpay);
	//EurPut_vega.Calc(para_volup);

	//std::vector<double> rs_vega = AmPut_vega.GetResult();
	//std::cout << "pv(AMERICAN)(for vega)" << rs_vega[0] << std::endl;
	//std::cout << "delta(for vega) " << rs_vega[1] << std::endl;
	//std::cout << "gamma(for vega) " << rs_vega[2] << std::endl;
	//std::cout << "spot(for vega) " << rs_vega[5] << std::endl;
	//std::cout << "BS vol(for vega)=" << volat_up.getBSVol(365 / 365, putpay.GetStrike()) << std::endl;
	//std::cout << "implied vol(for vega)=" << put_iv(spot, putstrike, rf_for_bs, div_for_bs, 0.1, 1.0, rs_vega[0]) << std::endl;
	////std::cout << "delta(for vega)(closed, Flatvol)=" << k_fmla_BSPutDelta(spot, putstrike, rf_for_bs, div_for_bs, volat_up.getBSVol(1.0, putpay.GetStrike()), 1.0) << std::endl;
	////std::cout << "gamma(for vega)(closed, Flatvol)=" << k_fmla_BSGamma(spot, putstrike, rf_for_bs, div_for_bs, volat_up.getBSVol(1.0, putpay.GetStrike()), 1.0) << std::endl;
	////std::cout << "\n";

	////std::cout << "vega(closed, Flatvol)=" << k_fmla_BSVega(spot, putstrike, rf_for_bs, div_for_bs, volat.getBSVol(1.0, putpay.GetStrike()), 1.0) / 100.0 << std::endl;
	////std::cout << "vega(for vega) " << rs_vega[0] - rs[0] << std::endl;


}



void test_asian(double spot, long numMC)
{

	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.54,1.64,1.72,1.77,1.8,1.85,1.88,1.90,1.92,1.92,1.92,1.92,1.92,1.92,1.92 };
	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.0191780821917808,0.038356164,0.057534247,0.084931507,0.104109589,0.123287671,0.142465753,0.167123288,0.252054795,0.334246575,0.419178082,0.504109589,0.580821918,0.665753425,0.747945205,0.832876712,0.915068493,1,1.084931507,1.167123288,1.252054795,1.334246575,1.419178082,1.504109589,1.583561644,1.668493151,1.750684932,1.835616438,1.917808219,	2.002739726	,2.087671233,2.169863014,2.254794521,2.336986301,2.421917808,2.506849315,2.583561644,2.668493151,2.750684932,2.835616438,2.917808219,3.002739726 };
	double arr_q[] = { 0.00000,0.00000,0.00000,2.57097, 2.09737,1.77112,1.53270, 1.30656 ,0.86631 ,5.16283, 4.11677,3.42318,3.34816,2.92103,2.60139,2.67302,2.43293,2.22630,2.26258, 2.10325,1.96057,3.00733,2.82806 ,2.66837, 2.53449,2.53675,2.41823 ,2.46233 ,2.35680 ,2.25685 ,2.27446 ,2.18831 ,2.10588 ,2.66917,2.57598,2.48871,2.41481,2.42319,2.35115 ,2.38170 ,2.31461 ,2.24914 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08493150,0.16712330,0.25205480,0.50410960,0.74794520,1.00000000,1.50410960,2.00273970,2.50684930,3.00273970 }; //10terms
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3095,0.2907,0.2719,0.2528,0.2331,0.2136,0.1816,0.1388,0.1081,0.105,0.1296,0.1606,0.1795,0.1892,0.1978,0.2057,0.2129 };
	double t1[] = { 0.2807,0.2619,0.2431,0.2239,0.2042,0.1845,0.1586,0.1285,0.1095,0.1027,0.1076,0.1306,0.1518,0.1634,0.172,0.1798,0.1869 };
	double t2[] = { 0.2621,0.2448,0.2274,0.2098,0.1916,0.172,0.1505,0.128,0.1131,0.1076,0.1125,0.1264,0.1413,0.1517,0.16,0.1676,0.1745 };
	double t3[] = { 0.2376,0.2219,0.2061,0.1902,0.1741,0.1576,0.141,0.1261,0.1169,0.1161,0.1215,0.1292,0.1371,0.1447,0.1517,0.1583,0.1644 };
	double t4[] = { 0.2286,0.214,0.1997,0.1854,0.1712,0.1571,0.1437,0.1319,0.1238,0.1211,0.1231,0.128,0.1339,0.14,0.146,0.1517,0.1571 };
	double t5[] = { 0.221,0.2077,0.1946,0.1819,0.1694,0.1574,0.1462,0.1366,0.1295,0.1261,0.1261,0.1287,0.1327,0.1374,0.1422,0.1469,0.1516 };
	double t6[] = { 0.2063,0.1949,0.1841,0.1737,0.1639,0.1549,0.1469,0.1403,0.1355,0.1329,0.1322,0.1331,0.1352,0.1379,0.1409,0.1442,0.1475 };
	double t7[] = { 0.1993,0.1897,0.1806,0.1721,0.1642,0.157,0.1506,0.1453,0.1411,0.1383,0.1368,0.1364,0.1369,0.1381,0.1398,0.1418,0.144 };
	double t8[] = { 0.1927,0.1844,0.1768,0.1696,0.1631,0.1572,0.1521,0.1477,0.1443,0.1418,0.1402,0.1394,0.1393,0.1398,0.1407,0.1419,0.1434 };
	double t9[] = { 0.1892,0.182,0.1754,0.1692,0.1637,0.1586,0.1542,0.1505,0.1474,0.145,0.1432,0.1421,0.1415,0.1414,0.1417,0.1423,0.1431 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);

	double refprice = 297.22;
	double rf_for_bs = r.getIntpRate(1.0);
	double div_for_bs = q.getIntpRate(1.0);
	double putstrike = 297.22*0.9;
	double callstrike = 297.22*1.1;
	double intpvol = volat.getInpVol(1.0, putstrike);

	MarketParam para(43340, spot, volat, r, q);

	PayoffPut putpay(putstrike);
	

	EuropeanOption EurPut(refprice, 43705, putpay);
	//vector<signed int> od{43432,43524,43613,43705};
	vector<signed int> od{ 43524,43705 };


	AsianOption Asian(refprice, 43705 ,putpay,od);
	Asian.Calc(para, numMC);
//	EurPut.Calc(para);

	/*EuropeanVanilla Vput(putpay,refprice, ed);

	EuropeanPut eput(refprice,putstrike,365);
	EuropeanCall ecall(refprice,callstrike,365);

	*/
	std::cout.precision(8);
	std::cout << std::fixed;

	//std::cout <<"eput:  s=" <<spot<<" "<< eput.pv(para)<<" vol " << volat.getInpVol(365/365,putstrike)<< "r " << r.getIntpRate(1.0) << " q " << q.getIntpRate(1.0) <<  std::endl;
	//std::cout <<"ecall: "<< ecall.pv(para)<< "s= " <<spot<<" "<<" vol " <<volat.getInpVol(365/365,callstrike)<< " r " << r.getIntpRate(1.0) << " q " << q.getIntpRate(1.0) << std::endl;
	std::vector<double> rs = Asian.GetResult();
	std::vector<double> euro_rs = EurPut.GetResult();

	std::cout << "   " <<"Asian " << "        "  << "european" << std::endl;
	std::cout << "pv " << rs[0] << "   " << euro_rs[0] <<std::endl;
	std::cout << "delta " << rs[1] << "  " << euro_rs[1]  << std::endl;
	std::cout << "gamma " << rs[2] <<  "  "  << euro_rs[2] << std::endl;
	std::cout << "spot " << rs[5] << std::endl;
	std::cout << "BS vol=" << volat.getBSVol(365 / 365, putpay.GetStrike()) << std::endl;
	std::cout << "implied vol=" << put_iv(spot, putstrike, rf_for_bs, div_for_bs, 0.1, 1.0, rs[0]) << std::endl;
	std::cout << "delta(closed, Flatvol)=" << k_fmla_BSPutDelta(spot, putstrike, rf_for_bs, div_for_bs, volat.getBSVol(1.0, putpay.GetStrike()), 1.0) << std::endl;
	std::cout << "gamma(closed, Flatvol)=" << k_fmla_BSGamma(spot, putstrike, rf_for_bs, div_for_bs, volat.getBSVol(1.0, putpay.GetStrike()), 1.0) << std::endl;

	std::cout << "\n";


}

void test_barrier_flat(double spot)
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.54,1.64,1.72,1.77,1.8,1.85,1.88,1.90,1.92,1.92,1.92,1.92,1.92,1.92,1.92 };
	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.0191780821917808,0.038356164,0.057534247,0.084931507,0.104109589,0.123287671,0.142465753,0.167123288,0.252054795,0.334246575,0.419178082,0.504109589,0.580821918,0.665753425,0.747945205,0.832876712,0.915068493,1,1.084931507,1.167123288,1.252054795,1.334246575,1.419178082,1.504109589,1.583561644,1.668493151,1.750684932,1.835616438,1.917808219,	2.002739726	,2.087671233,2.169863014,2.254794521,2.336986301,2.421917808,2.506849315,2.583561644,2.668493151,2.750684932,2.835616438,2.917808219,3.002739726 };
	double arr_q[] = { 0.00000,0.00000,0.00000,2.57097, 2.09737,1.77112,1.53270, 1.30656 ,0.86631 ,5.16283, 4.11677,3.42318,3.34816,2.92103,2.60139,2.67302,2.43293,2.22630,2.26258, 2.10325,1.96057,3.00733,2.82806 ,2.66837, 2.53449,2.53675,2.41823 ,2.46233 ,2.35680 ,2.25685 ,2.27446 ,2.18831 ,2.10588 ,2.66917,2.57598,2.48871,2.41481,2.42319,2.35115 ,2.38170 ,2.31461 ,2.24914 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08493150,0.16712330,0.25205480,0.50410960,0.74794520,1.00000000,1.50410960,2.00273970,2.50684930,3.00273970 }; //10terms
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3095,0.2907,0.2719,0.2528,0.2331,0.2136,0.1816,0.1388,0.1081,0.105,0.1296,0.1606,0.1795,0.1892,0.1978,0.2057,0.2129 };
	double t1[] = { 0.2807,0.2619,0.2431,0.2239,0.2042,0.1845,0.1586,0.1285,0.1095,0.1027,0.1076,0.1306,0.1518,0.1634,0.172,0.1798,0.1869 };
	double t2[] = { 0.2621,0.2448,0.2274,0.2098,0.1916,0.172,0.1505,0.128,0.1131,0.1076,0.1125,0.1264,0.1413,0.1517,0.16,0.1676,0.1745 };
	double t3[] = { 0.2376,0.2219,0.2061,0.1902,0.1741,0.1576,0.141,0.1261,0.1169,0.1161,0.1215,0.1292,0.1371,0.1447,0.1517,0.1583,0.1644 };
	double t4[] = { 0.2286,0.214,0.1997,0.1854,0.1712,0.1571,0.1437,0.1319,0.1238,0.1211,0.1231,0.128,0.1339,0.14,0.146,0.1517,0.1571 };
	double t5[] = { 0.221,0.2077,0.1946,0.1819,0.1694,0.1574,0.1462,0.1366,0.1295,0.1261,0.1261,0.1287,0.1327,0.1374,0.1422,0.1469,0.1516 };
	double t6[] = { 0.2063,0.1949,0.1841,0.1737,0.1639,0.1549,0.1469,0.1403,0.1355,0.1329,0.1322,0.1331,0.1352,0.1379,0.1409,0.1442,0.1475 };
	double t7[] = { 0.1993,0.1897,0.1806,0.1721,0.1642,0.157,0.1506,0.1453,0.1411,0.1383,0.1368,0.1364,0.1369,0.1381,0.1398,0.1418,0.144 };
	double t8[] = { 0.1927,0.1844,0.1768,0.1696,0.1631,0.1572,0.1521,0.1477,0.1443,0.1418,0.1402,0.1394,0.1393,0.1398,0.1407,0.1419,0.1434 };
	double t9[] = { 0.1892,0.182,0.1754,0.1692,0.1637,0.1586,0.1542,0.1505,0.1474,0.145,0.1432,0.1421,0.1415,0.1414,0.1417,0.1423,0.1431 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);

	double refprice = 297.22;
	double rf_for_bs = r.getIntpRate(1.0);
	double div_for_bs = q.getIntpRate(1.0);
	double putstrike = 297.22*0.9;
	double callstrike = 297.22*1.0;
	double upbarrier = 297.22*1.2;
	double intpvol = volat.getInpVol(1.0, putstrike);
	double rebate = 0.0;

	double const_vol = volat.getBSVol(1.0,297.22);
	double const_rf = r.getIntpRate(1.0);
	double const_q = q.getIntpRate(1.0);

	volat.set_const_vol(const_vol);
	r.set_const_rate(const_rf);
	q.set_const_rate(const_q);

	MarketParam para(0, spot, volat, r, q);

	PayoffBarCUO cuo(callstrike,upbarrier,rebate);
	
	BarrierOption Bar_cuo(refprice,365,cuo);

	Bar_cuo.Calc(para);

	//std::vector<double> rs = EurPut.GetResult();
	std::vector<double> rs_bar = Bar_cuo.GetResult();
	
	std::cout.precision(8);
	std::cout << std::fixed;

	std::cout << "pv " << rs_bar[0] << std::endl;
	std::cout << "delta " << rs_bar[1] << std::endl;
	std::cout << "gamma " << rs_bar[2] << std::endl;
	std::cout << "spot " << rs_bar[5] << std::endl;
	std::cout << "atm vol at ref " << volat.getBSVol(1.0, refprice) << std::endl;
	
	const double ds = 297.22*0.01;
	double clf_pv = k_fmla_Barrier(0, 0, 1, 297.22, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate);
	double clf_pv_up = k_fmla_Barrier(0, 0, 1, 297.22+ ds, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate);
	double clf_pv_down = k_fmla_Barrier(0, 0, 1, 297.22- ds, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate);

	std::cout << "closed form " << clf_pv << std::endl;
	std::cout << "closed form delta" << (clf_pv_up-clf_pv_down) /(ds*2.0) << std::endl;
	std::cout << "closed form gamma" << (clf_pv_up -2.0*clf_pv + clf_pv_down) /ds/ds << std::endl;

	std::cout << "\n";

	Vol volat_up(volat);
	volat_up.Vol_up(0.01);
	MarketParam para_volup(0, spot, volat_up, r, q);
	BarrierOption Bar_cuo_vega(refprice, 365, cuo);
	Bar_cuo_vega.Calc(para_volup);

	std::vector<double> rs_vega = Bar_cuo_vega.GetResult();
	std::cout << "\npv(for vega)" << rs_vega[0] << std::endl;
	std::cout << "delta(for vega) " << rs_vega[1] << std::endl;
	std::cout << "gamma(for vega) " << rs_vega[2] << std::endl;
	std::cout << "spot(for vega) " << rs_vega[5] << std::endl;
	std::cout << "BS vol(for vega)=" << volat_up.getBSVol(1.0, callstrike) << std::endl;
	std::cout << "\n";

	std::cout << "vega(closed, Flatvol)=" << k_fmla_Barrier(0, 0, 1, 297.22, callstrike, refprice*1.205355, const_rf, const_q, const_vol + 0.01, 1.0, rebate) - k_fmla_Barrier(0, 0, 1, 297.22, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate) << std::endl;
	std::cout << "vega(for vega) " << rs_vega[0] - rs_bar[0] << std::endl;


}

void test_barrier_LV(double spot)
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.54,1.64,1.72,1.77,1.8,1.85,1.88,1.90,1.92,1.92,1.92,1.92,1.92,1.92,1.92 };
	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.0191780821917808,0.038356164,0.057534247,0.084931507,0.104109589,0.123287671,0.142465753,0.167123288,0.252054795,0.334246575,0.419178082,0.504109589,0.580821918,0.665753425,0.747945205,0.832876712,0.915068493,1,1.084931507,1.167123288,1.252054795,1.334246575,1.419178082,1.504109589,1.583561644,1.668493151,1.750684932,1.835616438,1.917808219,	2.002739726	,2.087671233,2.169863014,2.254794521,2.336986301,2.421917808,2.506849315,2.583561644,2.668493151,2.750684932,2.835616438,2.917808219,3.002739726 };
	double arr_q[] = { 0.00000,0.00000,0.00000,2.57097, 2.09737,1.77112,1.53270, 1.30656 ,0.86631 ,5.16283, 4.11677,3.42318,3.34816,2.92103,2.60139,2.67302,2.43293,2.22630,2.26258, 2.10325,1.96057,3.00733,2.82806 ,2.66837, 2.53449,2.53675,2.41823 ,2.46233 ,2.35680 ,2.25685 ,2.27446 ,2.18831 ,2.10588 ,2.66917,2.57598,2.48871,2.41481,2.42319,2.35115 ,2.38170 ,2.31461 ,2.24914 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08493150,0.16712330,0.25205480,0.50410960,0.74794520,1.00000000,1.50410960,2.00273970,2.50684930,3.00273970 }; //10terms
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3095,0.2907,0.2719,0.2528,0.2331,0.2136,0.1816,0.1388,0.1081,0.105,0.1296,0.1606,0.1795,0.1892,0.1978,0.2057,0.2129 };
	double t1[] = { 0.2807,0.2619,0.2431,0.2239,0.2042,0.1845,0.1586,0.1285,0.1095,0.1027,0.1076,0.1306,0.1518,0.1634,0.172,0.1798,0.1869 };
	double t2[] = { 0.2621,0.2448,0.2274,0.2098,0.1916,0.172,0.1505,0.128,0.1131,0.1076,0.1125,0.1264,0.1413,0.1517,0.16,0.1676,0.1745 };
	double t3[] = { 0.2376,0.2219,0.2061,0.1902,0.1741,0.1576,0.141,0.1261,0.1169,0.1161,0.1215,0.1292,0.1371,0.1447,0.1517,0.1583,0.1644 };
	double t4[] = { 0.2286,0.214,0.1997,0.1854,0.1712,0.1571,0.1437,0.1319,0.1238,0.1211,0.1231,0.128,0.1339,0.14,0.146,0.1517,0.1571 };
	double t5[] = { 0.221,0.2077,0.1946,0.1819,0.1694,0.1574,0.1462,0.1366,0.1295,0.1261,0.1261,0.1287,0.1327,0.1374,0.1422,0.1469,0.1516 };
	double t6[] = { 0.2063,0.1949,0.1841,0.1737,0.1639,0.1549,0.1469,0.1403,0.1355,0.1329,0.1322,0.1331,0.1352,0.1379,0.1409,0.1442,0.1475 };
	double t7[] = { 0.1993,0.1897,0.1806,0.1721,0.1642,0.157,0.1506,0.1453,0.1411,0.1383,0.1368,0.1364,0.1369,0.1381,0.1398,0.1418,0.144 };
	double t8[] = { 0.1927,0.1844,0.1768,0.1696,0.1631,0.1572,0.1521,0.1477,0.1443,0.1418,0.1402,0.1394,0.1393,0.1398,0.1407,0.1419,0.1434 };
	double t9[] = { 0.1892,0.182,0.1754,0.1692,0.1637,0.1586,0.1542,0.1505,0.1474,0.145,0.1432,0.1421,0.1415,0.1414,0.1417,0.1423,0.1431 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);

	double refprice = 297.22;
	double rf_for_bs = r.getIntpRate(1.0);
	double div_for_bs = q.getIntpRate(1.0);
	double putstrike = 297.22*0.9;
	double callstrike = 297.22*1.0;
	double upbarrier = 297.22*1.2;
	double intpvol = volat.getInpVol(1.0, putstrike);
	double rebate = 0.0;

	double const_vol = volat.getBSVol(1.0, 297.22);
	double const_rf = r.getIntpRate(1.0);
	double const_q = q.getIntpRate(1.0);

	MarketParam para(0, spot, volat, r, q);

	PayoffBarCUO cuo(callstrike, upbarrier, rebate);

	BarrierOption Bar_cuo(refprice, 365, cuo);

	Bar_cuo.Calc(para);

	std::vector<double> rs_bar = Bar_cuo.GetResult();

	std::cout.precision(8);
	std::cout << std::fixed;

	std::cout << "\npv " << rs_bar[0] << std::endl;
	std::cout << "delta " << rs_bar[1] << std::endl;
	std::cout << "gamma " << rs_bar[2] << std::endl;
	std::cout << "spot " << rs_bar[5] << std::endl;
	std::cout << "atm vol at ref" << volat.getBSVol(1.0, refprice) << std::endl;

	const double ds = 297.22*0.01;
	double clf_pv = k_fmla_Barrier(0, 0, 1, 297.22, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate);
	double clf_pv_up = k_fmla_Barrier(0, 0, 1, 297.22 + ds, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate);
	double clf_pv_down = k_fmla_Barrier(0, 0, 1, 297.22 - ds, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate);

	std::cout << "closed form " << clf_pv << std::endl;
	std::cout << "closed form delta" << (clf_pv_up - clf_pv_down) / (ds*2.0) << std::endl;
	std::cout << "closed form gamma" << (clf_pv_up - 2.0*clf_pv + clf_pv_down) / ds / ds << std::endl;
	std::cout << "\n";

	Vol volat_up(volat);
	volat_up.Vol_up(0.01);
	MarketParam para_volup(0, spot, volat_up, r, q);
	BarrierOption Bar_cuo_vega(refprice, 365, cuo);
	Bar_cuo_vega.Calc(para_volup);

	std::vector<double> rs_vega = Bar_cuo_vega.GetResult();
	std::cout << "pv(for vega)" << rs_vega[0] << std::endl;
	std::cout << "delta(for vega) " << rs_vega[1] << std::endl;
	std::cout << "gamma(for vega) " << rs_vega[2] << std::endl;
	std::cout << "spot(for vega) " << rs_vega[5] << std::endl;
	std::cout << "BS vol(for vega)=" << volat_up.getBSVol(1.0, callstrike) << std::endl;
	std::cout << "\n";

	std::cout << "vega(closed, Flatvol)=" << k_fmla_Barrier(0, 0, 1, 297.22, callstrike, refprice*1.205355, const_rf, const_q, const_vol+0.01, 1.0, rebate)- k_fmla_Barrier(0, 0, 1, 297.22, callstrike, refprice*1.205355, const_rf, const_q, const_vol, 1.0, rebate) << std::endl;
	std::cout << "vega(for vega) " << rs_vega[0] - rs_bar[0] << std::endl;

}


MarketParameters set_paras()
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.79000,1.81000,1.74000,1.69000,1.65000,1.59000,1.55000,1.52000,1.50000,1.50000,1.50000,1.50000,1.50000,1.50000 ,1.50000 };

	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.01918000 ,	0.03836000 ,	0.05753000 ,	0.07671000 	,0.09589000 ,	0.11507000 	,0.13425000 ,	0.16438000 ,	0.24658000, 	0.32877000 ,	0.41096000 	,0.49315000 	,0.57534000 ,	0.65753000 ,	0.73973000 ,	0.82192000 ,	0.90411000 ,	0.98630000 ,	1.06849000 ,	1.15068000 	,1.23288000 ,	1.31507000 ,	1.39726000 ,	1.47945000 	,1.56164000 	,1.64384000 ,	1.72603000 ,	1.80822000 ,	1.89041000 ,	1.97260000 ,	2.05479000 ,2.13699000 	,2.21918000, 	2.30137000 ,	2.38356000 	,2.46575000 ,	2.54795000 ,	2.63014000 ,	2.71233000 ,	2.79452000 ,	2.87671000, 	2.95890000 };
	double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,2.33330,1.94438,1.666590,1.361110,0.907370,0.680540,4.321890,3.601590,3.087090,3.035340,2.698050,2.428250,2.517120,2.307370,2.129880,2.172190,2.027360,1.900660,2.957910,2.793580,2.646550,2.648480,2.522360,2.407710,2.456450,2.354100,2.259940,2.278180,2.193800,2.115460,2.744900,2.653400,2.567800,2.571480,2.493550,2.420220,2.451880,2.383780 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08333,0.16667,0.250000,0.500000,0.750000,1.000000,1.500000,2.000000,2.500000,3.000000 };
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3484,	0.3310,	0.3135,	0.2960,	0.2778,	0.2577, 0.2279, 0.1924, 0.1641, 0.1532, 0.1618, 0.1729, 0.1807, 0.1883, 0.1968,	0.2034,	0.2090 };
	double t1[] = { 0.3182,	0.3002,	0.2822,	0.2640,	0.2454,	0.2252,	0.2010,	0.1770,	0.1568,	0.1452,	0.1429,	0.1484,	0.1565,	0.1645,	0.1722,	0.1792,	0.1855 };
	double t2[] = { 0.3069,	0.2883,	0.2697,	0.2509,	0.2318,	0.2118,	0.1907,	0.1699,	0.1527,	0.1432,	0.1423,	0.1476,	0.1554,	0.1634,	0.1711,	0.1782,	0.1849 };
	double t3[] = { 0.2709,	0.2543,	0.2380,	0.2218,	0.2058,	0.1900,	0.1748,	0.1607,	0.1491,	0.1418,	0.1395,	0.1411,	0.1451,	0.1502,	0.1555,	0.1609,	0.1662 };
	double t4[] = { 0.2441,	0.2300,	0.2163,	0.2029,	0.1901,	0.1778,	0.1664,	0.1562,	0.1479,	0.1423,	0.1396,	0.1394,	0.1410,	0.1437,	0.1470,	0.1506,	0.1543 };
	double t5[] = { 0.2277,	0.2156,	0.2041,	0.1930,	0.1826,	0.1727,	0.1636,	0.1556,	0.1490,	0.1440,	0.1409,	0.1395,	0.1396,	0.1406,	0.1424,	0.1446,	0.1471 };
	double t6[] = { 0.2073, 0.1981,	0.1895,	0.1814,	0.1739,	0.1669,	0.1606,	0.1550,	0.1503,	0.1465,	0.1436,	0.1416,	0.1404,	0.1400,	0.1401,	0.1407,	0.1417 };
	double t7[] = { 0.1944,	0.1871,	0.1802,	0.1739,	0.1681,	0.1629,	0.1581,	0.1539,	0.1503,	0.1473,	0.1449,	0.1431,	0.1418,	0.1410,	0.1406,	0.1405,	0.1407 };
	double t8[] = { 0.1876,	0.1815,	0.1760,	0.1709,	0.1662,	0.1619,	0.1580,	0.1546,	0.1516,	0.1490,	0.1468,	0.1451,	0.1436,	0.1426,	0.1418,	0.1413,	0.1410 };
	double t9[] = { 0.1821,	0.1770,	0.1723,	0.1680,	0.1641,	0.1606,	0.1573,	0.1545,	0.1520,	0.1497,	0.1478,	0.1462,	0.1449,	0.1438,	0.1429,	0.1422,	0.1418 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);

	//***상품타입변경
	double spot = 297.22;
	signed int vd = 43340;
	
	//return MarketParameters();
	return MarketParameters(vd, spot, volat, r, q);
}

MarketParameters set_paras_spot1()
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.79000,1.81000,1.74000,1.69000,1.65000,1.59000,1.55000,1.52000,1.50000,1.50000,1.50000,1.50000,1.50000,1.50000 ,1.50000 };

	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.01918000 ,	0.03836000 ,	0.05753000 ,	0.07671000 	,0.09589000 ,	0.11507000 	,0.13425000 ,	0.16438000 ,	0.24658000, 	0.32877000 ,	0.41096000 	,0.49315000 	,0.57534000 ,	0.65753000 ,	0.73973000 ,	0.82192000 ,	0.90411000 ,	0.98630000 ,	1.06849000 ,	1.15068000 	,1.23288000 ,	1.31507000 ,	1.39726000 ,	1.47945000 	,1.56164000 	,1.64384000 ,	1.72603000 ,	1.80822000 ,	1.89041000 ,	1.97260000 ,	2.05479000 ,2.13699000 	,2.21918000, 	2.30137000 ,	2.38356000 	,2.46575000 ,	2.54795000 ,	2.63014000 ,	2.71233000 ,	2.79452000 ,	2.87671000, 	2.95890000 };
	double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,2.33330,1.94438,1.666590,1.361110,0.907370,0.680540,4.321890,3.601590,3.087090,3.035340,2.698050,2.428250,2.517120,2.307370,2.129880,2.172190,2.027360,1.900660,2.957910,2.793580,2.646550,2.648480,2.522360,2.407710,2.456450,2.354100,2.259940,2.278180,2.193800,2.115460,2.744900,2.653400,2.567800,2.571480,2.493550,2.420220,2.451880,2.383780 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08333,0.16667,0.250000,0.500000,0.750000,1.000000,1.500000,2.000000,2.500000,3.000000 };
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3484,	0.3310,	0.3135,	0.2960,	0.2778,	0.2577, 0.2279, 0.1924, 0.1641, 0.1532, 0.1618, 0.1729, 0.1807, 0.1883, 0.1968,	0.2034,	0.2090 };
	double t1[] = { 0.3182,	0.3002,	0.2822,	0.2640,	0.2454,	0.2252,	0.2010,	0.1770,	0.1568,	0.1452,	0.1429,	0.1484,	0.1565,	0.1645,	0.1722,	0.1792,	0.1855 };
	double t2[] = { 0.3069,	0.2883,	0.2697,	0.2509,	0.2318,	0.2118,	0.1907,	0.1699,	0.1527,	0.1432,	0.1423,	0.1476,	0.1554,	0.1634,	0.1711,	0.1782,	0.1849 };
	double t3[] = { 0.2709,	0.2543,	0.2380,	0.2218,	0.2058,	0.1900,	0.1748,	0.1607,	0.1491,	0.1418,	0.1395,	0.1411,	0.1451,	0.1502,	0.1555,	0.1609,	0.1662 };
	double t4[] = { 0.2441,	0.2300,	0.2163,	0.2029,	0.1901,	0.1778,	0.1664,	0.1562,	0.1479,	0.1423,	0.1396,	0.1394,	0.1410,	0.1437,	0.1470,	0.1506,	0.1543 };
	double t5[] = { 0.2277,	0.2156,	0.2041,	0.1930,	0.1826,	0.1727,	0.1636,	0.1556,	0.1490,	0.1440,	0.1409,	0.1395,	0.1396,	0.1406,	0.1424,	0.1446,	0.1471 };
	double t6[] = { 0.2073, 0.1981,	0.1895,	0.1814,	0.1739,	0.1669,	0.1606,	0.1550,	0.1503,	0.1465,	0.1436,	0.1416,	0.1404,	0.1400,	0.1401,	0.1407,	0.1417 };
	double t7[] = { 0.1944,	0.1871,	0.1802,	0.1739,	0.1681,	0.1629,	0.1581,	0.1539,	0.1503,	0.1473,	0.1449,	0.1431,	0.1418,	0.1410,	0.1406,	0.1405,	0.1407 };
	double t8[] = { 0.1876,	0.1815,	0.1760,	0.1709,	0.1662,	0.1619,	0.1580,	0.1546,	0.1516,	0.1490,	0.1468,	0.1451,	0.1436,	0.1426,	0.1418,	0.1413,	0.1410 };
	double t9[] = { 0.1821,	0.1770,	0.1723,	0.1680,	0.1641,	0.1606,	0.1573,	0.1545,	0.1520,	0.1497,	0.1478,	0.1462,	0.1449,	0.1438,	0.1429,	0.1422,	0.1418 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);

	double spot = 297.22-1.0;

	signed int vd = 43340;

	return MarketParameters(vd, spot, volat, r, q);
}

MarketParameters set_paras_flat()
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.79000,1.81000,1.74000,1.69000,1.65000,1.59000,1.55000,1.52000,1.50000,1.50000,1.50000,1.50000,1.50000,1.50000 ,1.50000 };
	//double arr_r[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,0 };

	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.01918000 ,	0.03836000 ,	0.05753000 ,	0.07671000 	,0.09589000 ,	0.11507000 	,0.13425000 ,	0.16438000 ,	0.24658000, 	0.32877000 ,	0.41096000 	,0.49315000 	,0.57534000 ,	0.65753000 ,	0.73973000 ,	0.82192000 ,	0.90411000 ,	0.98630000 ,	1.06849000 ,	1.15068000 	,1.23288000 ,	1.31507000 ,	1.39726000 ,	1.47945000 	,1.56164000 	,1.64384000 ,	1.72603000 ,	1.80822000 ,	1.89041000 ,	1.97260000 ,	2.05479000 ,2.13699000 	,2.21918000, 	2.30137000 ,	2.38356000 	,2.46575000 ,	2.54795000 ,	2.63014000 ,	2.71233000 ,	2.79452000 ,	2.87671000, 	2.95890000 };
	double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,2.33330,1.94438,1.666590,1.361110,0.907370,0.680540,4.321890,3.601590,3.087090,3.035340,2.698050,2.428250,2.517120,2.307370,2.129880,2.172190,2.027360,1.900660,2.957910,2.793580,2.646550,2.648480,2.522360,2.407710,2.456450,2.354100,2.259940,2.278180,2.193800,2.115460,2.744900,2.653400,2.567800,2.571480,2.493550,2.420220,2.451880,2.383780 };
	//double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08333,0.16667,0.250000,0.500000,0.750000,1.000000,1.500000,2.000000,2.500000,3.000000 };
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	//double t0[] = { 0.3484,	0.3310,	0.3135,	0.2960,	0.2778,	0.2577, 0.2279, 0.1924, 0.1641, 0.1532, 0.1618, 0.1729, 0.1807, 0.1883, 0.1968,	0.2034,	0.2090 };
	//double t1[] = { 0.3182,	0.3002,	0.2822,	0.2640,	0.2454,	0.2252,	0.2010,	0.1770,	0.1568,	0.1452,	0.1429,	0.1484,	0.1565,	0.1645,	0.1722,	0.1792,	0.1855 };
	//double t2[] = { 0.3069,	0.2883,	0.2697,	0.2509,	0.2318,	0.2118,	0.1907,	0.1699,	0.1527,	0.1432,	0.1423,	0.1476,	0.1554,	0.1634,	0.1711,	0.1782,	0.1849 };
	//double t3[] = { 0.2709,	0.2543,	0.2380,	0.2218,	0.2058,	0.1900,	0.1748,	0.1607,	0.1491,	0.1418,	0.1395,	0.1411,	0.1451,	0.1502,	0.1555,	0.1609,	0.1662 };
	//double t4[] = { 0.2441,	0.2300,	0.2163,	0.2029,	0.1901,	0.1778,	0.1664,	0.1562,	0.1479,	0.1423,	0.1396,	0.1394,	0.1410,	0.1437,	0.1470,	0.1506,	0.1543 };
	//double t5[] = { 0.2277,	0.2156,	0.2041,	0.1930,	0.1826,	0.1727,	0.1636,	0.1556,	0.1490,	0.1440,	0.1409,	0.1395,	0.1396,	0.1406,	0.1424,	0.1446,	0.1471 };
	//double t6[] = { 0.2073, 0.1981,	0.1895,	0.1814,	0.1739,	0.1669,	0.1606,	0.1550,	0.1503,	0.1465,	0.1436,	0.1416,	0.1404,	0.1400,	0.1401,	0.1407,	0.1417 };
	//double t7[] = { 0.1944,	0.1871,	0.1802,	0.1739,	0.1681,	0.1629,	0.1581,	0.1539,	0.1503,	0.1473,	0.1449,	0.1431,	0.1418,	0.1410,	0.1406,	0.1405,	0.1407 };
	//double t8[] = { 0.1876,	0.1815,	0.1760,	0.1709,	0.1662,	0.1619,	0.1580,	0.1546,	0.1516,	0.1490,	0.1468,	0.1451,	0.1436,	0.1426,	0.1418,	0.1413,	0.1410 };
	//double t9[] = { 0.1821,	0.1770,	0.1723,	0.1680,	0.1641,	0.1606,	0.1573,	0.1545,	0.1520,	0.1497,	0.1478,	0.1462,	0.1449,	0.1438,	0.1429,	0.1422,	0.1418 };

	double t0[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t1[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t2[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t3[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t4[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t5[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t6[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t7[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t8[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };
	double t9[] = { 0.1490,	0.1490,	0.1490,	0.1490,	0.1490,	0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490, 0.1490,	0.1490,	0.1490 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);


	//double intpvol = volat.getInpVol(1.0, putstrike);  //flat interpolation use .getBSVol

	//상품타입변경
	//bool is_flat_vol = true;
	//int hitflag = 1;
	//long nM = 30000;

	//***상품타입변경
	double spot = 297.22;
	signed int vd = 43340;

	return MarketParameters(vd, spot, volat, r, q);
}

MarketParam set_para_flat()
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.79000,1.81000,1.74000,1.69000,1.65000,1.59000,1.55000,1.52000,1.50000,1.50000,1.50000,1.50000,1.50000,1.50000 ,1.50000 };
	//double arr_r[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,0 };

	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.01918000 ,	0.03836000 ,	0.05753000 ,	0.07671000 	,0.09589000 ,	0.11507000 	,0.13425000 ,	0.16438000 ,	0.24658000, 	0.32877000 ,	0.41096000 	,0.49315000 	,0.57534000 ,	0.65753000 ,	0.73973000 ,	0.82192000 ,	0.90411000 ,	0.98630000 ,	1.06849000 ,	1.15068000 	,1.23288000 ,	1.31507000 ,	1.39726000 ,	1.47945000 	,1.56164000 	,1.64384000 ,	1.72603000 ,	1.80822000 ,	1.89041000 ,	1.97260000 ,	2.05479000 ,2.13699000 	,2.21918000, 	2.30137000 ,	2.38356000 	,2.46575000 ,	2.54795000 ,	2.63014000 ,	2.71233000 ,	2.79452000 ,	2.87671000, 	2.95890000 };
	//double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,2.33330,1.94438,1.666590,1.361110,0.907370,0.680540,4.321890,3.601590,3.087090,3.035340,2.698050,2.428250,2.517120,2.307370,2.129880,2.172190,2.027360,1.900660,2.957910,2.793580,2.646550,2.648480,2.522360,2.407710,2.456450,2.354100,2.259940,2.278180,2.193800,2.115460,2.744900,2.653400,2.567800,2.571480,2.493550,2.420220,2.451880,2.383780 };
	double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08333,0.16667,0.250000,0.500000,0.750000,1.000000,1.500000,2.000000,2.500000,3.000000 };
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3484,	0.3310,	0.3135,	0.2960,	0.2778,	0.2577, 0.2279, 0.1924, 0.1641, 0.1532, 0.1618, 0.1729, 0.1807, 0.1883, 0.1968,	0.2034,	0.2090 };
	double t1[] = { 0.3182,	0.3002,	0.2822,	0.2640,	0.2454,	0.2252,	0.2010,	0.1770,	0.1568,	0.1452,	0.1429,	0.1484,	0.1565,	0.1645,	0.1722,	0.1792,	0.1855 };
	double t2[] = { 0.3069,	0.2883,	0.2697,	0.2509,	0.2318,	0.2118,	0.1907,	0.1699,	0.1527,	0.1432,	0.1423,	0.1476,	0.1554,	0.1634,	0.1711,	0.1782,	0.1849 };
	double t3[] = { 0.2709,	0.2543,	0.2380,	0.2218,	0.2058,	0.1900,	0.1748,	0.1607,	0.1491,	0.1418,	0.1395,	0.1411,	0.1451,	0.1502,	0.1555,	0.1609,	0.1662 };
	double t4[] = { 0.2441,	0.2300,	0.2163,	0.2029,	0.1901,	0.1778,	0.1664,	0.1562,	0.1479,	0.1423,	0.1396,	0.1394,	0.1410,	0.1437,	0.1470,	0.1506,	0.1543 };
	double t5[] = { 0.2277,	0.2156,	0.2041,	0.1930,	0.1826,	0.1727,	0.1636,	0.1556,	0.1490,	0.1440,	0.1409,	0.1395,	0.1396,	0.1406,	0.1424,	0.1446,	0.1471 };
	double t6[] = { 0.2073, 0.1981,	0.1895,	0.1814,	0.1739,	0.1669,	0.1606,	0.1550,	0.1503,	0.1465,	0.1436,	0.1416,	0.1404,	0.1400,	0.1401,	0.1407,	0.1417 };
	double t7[] = { 0.1944,	0.1871,	0.1802,	0.1739,	0.1681,	0.1629,	0.1581,	0.1539,	0.1503,	0.1473,	0.1449,	0.1431,	0.1418,	0.1410,	0.1406,	0.1405,	0.1407 };
	double t8[] = { 0.1876,	0.1815,	0.1760,	0.1709,	0.1662,	0.1619,	0.1580,	0.1546,	0.1516,	0.1490,	0.1468,	0.1451,	0.1436,	0.1426,	0.1418,	0.1413,	0.1410 };
	double t9[] = { 0.1821,	0.1770,	0.1723,	0.1680,	0.1641,	0.1606,	0.1573,	0.1545,	0.1520,	0.1497,	0.1478,	0.1462,	0.1449,	0.1438,	0.1429,	0.1422,	0.1418 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);


	//double intpvol = volat.getInpVol(1.0, putstrike);  //flat interpolation use .getBSVol

	//상품타입변경
	//bool is_flat_vol = true;
	//int hitflag = 1;
	//long nM = 30000;

	//***상품타입변경
	double spot = 297.22;
	signed int vd = 43340;

	return MarketParam(vd, spot, volat, r, q);
}

MarketParam set_para()
{
	double arr_rts[] = { 0.002739726,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,7,10,15,20 };
	double arr_r[] = { 1.79000,1.81000,1.74000,1.69000,1.65000,1.59000,1.55000,1.52000,1.50000,1.50000,1.50000,1.50000,1.50000,1.50000 ,1.50000 };

	for (int i = 0; i<15; i++)
		arr_r[i] /= 100.0;

	double arr_qts[] = { 0.01918000 ,	0.03836000 ,	0.05753000 ,	0.07671000 	,0.09589000 ,	0.11507000 	,0.13425000 ,	0.16438000 ,	0.24658000, 	0.32877000 ,	0.41096000 	,0.49315000 	,0.57534000 ,	0.65753000 ,	0.73973000 ,	0.82192000 ,	0.90411000 ,	0.98630000 ,	1.06849000 ,	1.15068000 	,1.23288000 ,	1.31507000 ,	1.39726000 ,	1.47945000 	,1.56164000 	,1.64384000 ,	1.72603000 ,	1.80822000 ,	1.89041000 ,	1.97260000 ,	2.05479000 ,2.13699000 	,2.21918000, 	2.30137000 ,	2.38356000 	,2.46575000 ,	2.54795000 ,	2.63014000 ,	2.71233000 ,	2.79452000 ,	2.87671000, 	2.95890000 };
	double arr_q[] = { 0.00000,0.00000,0.00000,0.00000,2.33330,1.94438,1.666590,1.361110,0.907370,0.680540,4.321890,3.601590,3.087090,3.035340,2.698050,2.428250,2.517120,2.307370,2.129880,2.172190,2.027360,1.900660,2.957910,2.793580,2.646550,2.648480,2.522360,2.407710,2.456450,2.354100,2.259940,2.278180,2.193800,2.115460,2.744900,2.653400,2.567800,2.571480,2.493550,2.420220,2.451880,2.383780 };
	for (int i = 0; i<42; i++)
		arr_q[i] /= 100.0;

	//double spot=297.22;
	double vol_term[] = { 0.08333,0.16667,0.250000,0.500000,0.750000,1.000000,1.500000,2.000000,2.500000,3.000000 };
	double vol_strike[] = { 178.332,193.193,208.054,222.915,237.776,252.637,267.498,282.359,297.22,312.081,326.942,341.803,356.664,371.525,386.386,401.247,416.108 };//17spot
	Vol volat(10, 17);
	double t0[] = { 0.3484,	0.3310,	0.3135,	0.2960,	0.2778,	0.2577, 0.2279, 0.1924, 0.1641, 0.1532, 0.1618, 0.1729, 0.1807, 0.1883, 0.1968,	0.2034,	0.2090 };
	double t1[] = { 0.3182,	0.3002,	0.2822,	0.2640,	0.2454,	0.2252,	0.2010,	0.1770,	0.1568,	0.1452,	0.1429,	0.1484,	0.1565,	0.1645,	0.1722,	0.1792,	0.1855 };
	double t2[] = { 0.3069,	0.2883,	0.2697,	0.2509,	0.2318,	0.2118,	0.1907,	0.1699,	0.1527,	0.1432,	0.1423,	0.1476,	0.1554,	0.1634,	0.1711,	0.1782,	0.1849 };
	double t3[] = { 0.2709,	0.2543,	0.2380,	0.2218,	0.2058,	0.1900,	0.1748,	0.1607,	0.1491,	0.1418,	0.1395,	0.1411,	0.1451,	0.1502,	0.1555,	0.1609,	0.1662 };
	double t4[] = { 0.2441,	0.2300,	0.2163,	0.2029,	0.1901,	0.1778,	0.1664,	0.1562,	0.1479,	0.1423,	0.1396,	0.1394,	0.1410,	0.1437,	0.1470,	0.1506,	0.1543 };
	double t5[] = { 0.2277,	0.2156,	0.2041,	0.1930,	0.1826,	0.1727,	0.1636,	0.1556,	0.1490,	0.1440,	0.1409,	0.1395,	0.1396,	0.1406,	0.1424,	0.1446,	0.1471 };
	double t6[] = { 0.2073, 0.1981,	0.1895,	0.1814,	0.1739,	0.1669,	0.1606,	0.1550,	0.1503,	0.1465,	0.1436,	0.1416,	0.1404,	0.1400,	0.1401,	0.1407,	0.1417 };
	double t7[] = { 0.1944,	0.1871,	0.1802,	0.1739,	0.1681,	0.1629,	0.1581,	0.1539,	0.1503,	0.1473,	0.1449,	0.1431,	0.1418,	0.1410,	0.1406,	0.1405,	0.1407 };
	double t8[] = { 0.1876,	0.1815,	0.1760,	0.1709,	0.1662,	0.1619,	0.1580,	0.1546,	0.1516,	0.1490,	0.1468,	0.1451,	0.1436,	0.1426,	0.1418,	0.1413,	0.1410 };
	double t9[] = { 0.1821,	0.1770,	0.1723,	0.1680,	0.1641,	0.1606,	0.1573,	0.1545,	0.1520,	0.1497,	0.1478,	0.1462,	0.1449,	0.1438,	0.1429,	0.1422,	0.1418 };

	volat.set_vol_strike(vol_strike, 17);
	volat.set_vol_term(vol_term, 10);

	volat.set_volsurface_by_term(t0, 0, 17);
	volat.set_volsurface_by_term(t1, 1, 17);
	volat.set_volsurface_by_term(t2, 2, 17);
	volat.set_volsurface_by_term(t3, 3, 17);
	volat.set_volsurface_by_term(t4, 4, 17);
	volat.set_volsurface_by_term(t5, 5, 17);
	volat.set_volsurface_by_term(t6, 6, 17);
	volat.set_volsurface_by_term(t7, 7, 17);
	volat.set_volsurface_by_term(t8, 8, 17);
	volat.set_volsurface_by_term(t9, 9, 17);

	std::vector<double> v_rts(std::begin(arr_rts), std::end(arr_rts));
	std::vector<double> v_r(std::begin(arr_r), std::end(arr_r));
	std::vector<double> v_qts(std::begin(arr_qts), std::end(arr_qts));
	std::vector<double> v_q(std::begin(arr_q), std::end(arr_q));

	Rate r(v_r, v_rts);
	Rate q(v_q, v_qts);


	//double intpvol = volat.getInpVol(1.0, putstrike);  //flat interpolation use .getBSVol

	//상품타입변경
	//bool is_flat_vol = true;
	//int hitflag = 1;
	//long nM = 30000;

	//***상품타입변경
	double spot = 297.22;
	signed int vd = 43340;


	return MarketParam(vd, spot, volat, r, q);
}



//int getDelta(double target, vector<double>& old, unsigned int& init_i)
//{
//	if (target <= old.front())
//		return (init_i = 0);
//	if (target >= old.back())
//		return (init_i = old.size()-1);
//
//	if (old[init_i] <= target && target <old[init_i + 1]) {
//		if (target - old[init_i] <old[init_i + 1] - target) {
//			return init_i;
//		}
//		else {
//			return (init_i = init_i + 1);
//		}
//	}
//
//	int i = init_i;
//	while (1) {   //향후 이부분 개선 
//		i = i + 1;
//		if (i < old.size() - 1) {
//			if (old[i] <= target && target < old[i + 1]) {
//				if (target - old[i] < old[i + 1] - target) {
//					return (init_i = i);
//				}
//				else {
//					return (init_i = i + 1);
//				}
//			}
//		}
//
//		int j = init_i - (i - init_i);
//		if (j >= 0) {
//			if (old[j] <= target && target <old[j + 1]) {
//				if (target - old[j] <old[j + 1] - target) {
//					return (init_i = j);
//				}
//				else {
//					return (init_i = j + 1);
//				}
//			}
//		}
//	}
//
//	throw std::logic_error("find_index_spot2 - interpolaton fail :findnearestindex, vol strike");
//	return -1;
//
//}
//
//unsigned int getIndex(double target, double* px, int i_min, int i_max, unsigned int init_i)
//{
//	if (target <= px[0])
//		return (init_i = 0);
//	if (target >= px[i_max])
//		return (init_i = i_max);
//
//	if (px[init_i] <= target && target <px[init_i + 1]) {
//		if (target - px[init_i] <px[init_i + 1] - target) {
//			return init_i;
//		}
//		else {
//			return (init_i = init_i + 1);
//		}
//	}
//
//	int i = init_i;
//	while (1) {   //향후 이부분 개선 
//		i = i + 1;
//		if (i < i_max) {
//			if (px[i] <= target && target < px[i + 1]) {
//				if (target - px[i] < px[i + 1] - target) {
//					return (init_i = i);
//				}
//				else {
//					return (init_i = i + 1);
//				}
//			}
//		}
//
//		int j = init_i - (i - init_i);
//		if (j >= 0) {
//			if (px[j] <= target && target <px[j + 1]) {
//				if (target - px[j] < px[j + 1] - target) {
//					return (init_i = j);
//				}
//				else {
//					return (init_i = j + 1);
//				}
//			}
//		}
//	}
//
//	throw std::logic_error("find_index_spot2 - interpolaton fail :findnearestindex, vol strike");
//	return -1;
//}


vector<double> get_a_path_from_csv(char* s)
{
	ifstream ifs;
	string str_buf;
	vector<double> apath;
	ifs.open(s);
	while (!ifs.eof()) {
		getline(ifs, str_buf, '\n');
		apath.push_back(stod(str_buf));
	}
	ifs.close();
	return apath;
}

void test_vanilla_iofile()
{
	MarketParameters paras_file = init_paras_file("vol20200129new.csv", "rate20200129new.csv", "div20200129new.csv");
	//discrete dividend
	MarketParameters paras_file_discrete = init_paras_file("vol20200129new.csv", "rate20200129new.csv", "div20200129new.csv", "divDiscrete20200129.csv");
	
	signed int exd6M = paras_file.get_vdate() + 180; //6month
	double refprice = paras_file.get_spot();
	double putstrike90 = refprice*0.9;
	double putstrike = refprice;

	PayoffPut putpay90(putstrike90);
	PayoffPut putpay(putstrike);

	EuropeanOption EurPut6M90(refprice, exd6M, putpay90);
	EuropeanOption EurPut6M100(refprice, exd6M, putpay);

	EurPut6M90.CalcBS(paras_file);
	cout << "\nEurPut1YATM.CalcBS(paras_file);" << endl;
	EurPut6M90.PrintResult();

	EurPut6M90.Calc(paras_file); 
	cout << "\nEurPut1YATM.Calc(paras_file);" << endl;
	EurPut6M90.PrintResult();
	paras_file.print();

	//test fdm discrete dividend
	EurPut6M90.Calc_discrete(paras_file_discrete);
	cout << "\nEurPut6M90.Calc_discrete(paras_file_discrete)" << endl;
	EurPut6M90.PrintResult();
	paras_file.print();

	//test mc continuous dividend
	EurPut6M90.CalcMC(paras_file, 50000);
	cout << "\nEurPut6M90.CalcMC(paras_file, 20000);" << endl;
	EurPut6M90.PrintResult();

	//test mc discrete dividend
	EurPut6M90.CalcMC_discrete(paras_file_discrete,50000);
	cout << "\nEurPut6M90.CalcMC_discrete(paras_file_discrete,20000);" << endl;
	EurPut6M90.PrintResult();
	paras_file.print();

	//vector<double> apath;
	//apath = get_a_path_from_csv("14th.csv");
	//EurPut6M90.Simulation2(paras_file,10, true);
	//EurPut6M100.Simulation2(paras_file, 10, true);

}

void test_vanilla_final()
{
	MarketParameters paras;
	MarketParameters paras_volup;
	paras = set_paras();
	paras_volup = set_paras();

	signed int exd = paras.get_vdate() + 365; //1y vanilla
	signed int exd6m = paras.get_vdate() + 180; //1y vanilla

	double refprice = 297.22;
	double putstrike = refprice;
	double putstrike90 = refprice*0.9;
	PayoffPut putpay(putstrike);
	PayoffPut put90(putstrike90);
	
	EuropeanOption EurPut1YATM(refprice, exd, putpay);
	EuropeanOption EurPut1Y90(refprice, exd, put90);
	EuropeanOption EurPut6MATM(refprice, exd6m, putpay);

	//EurPut1YATM.Calc(paras);
	//cout << "\nEurPut1YATM.Calc(paras);" << endl;
	//EurPut1YATM.PrintResult();

	//EurPut1YATM.Calc_discrete(paras);
	//cout << "\n	EurPut1YATM.Calc_discrete(paras);" << endl;
	//EurPut1YATM.PrintResult();

	//EurPut1YATM.CalcMC(paras,100000);
	//cout << "\nEurPut1YATM.CalcMC(paras);" << endl;
	//EurPut1YATM.PrintResult();

	//double bsval=EurPut1YATM.CalcBS(paras);
	//cout << "\nEurPut1YATM.CalcBS(paras);" << endl;
	//EurPut1YATM.PrintResult();
	//cout << "BS implied vol= " << EurPut1YATM.BSiv(paras, bsval);

}


void test_autocall_final_iofile()
{
	MarketParameters paras_file = init_paras_file("vol20200129new.csv", "rate20200129new.csv", "div20200129new.csv");
	MarketParameters paras_file_volup = init_paras_file("vol20200129new.csv", "rate20200129new.csv", "div20200129new.csv");
	MarketParameters paras_20180828 = init_paras_file("vol20180828.csv", "rate20180828.csv", "div20180828.csv");
	MarketParameters paras_file_div = init_paras_file("vol20200129new.csv", "rate20200129new.csv", "div20200129new.csv", "divDiscrete20200129.csv");

	AutocallOption prot("autocall20200129protected.csv");
	AutocallOption notprot("autocall20200129notprotected.csv");
	//AutocallOption notprot20180828("autocall20180828.csv");


	paras_file_volup.reset_Ivol_up();
	//paras_file.print();
	//cout << "volup\n";
	//paras_file_volup.print();

	//ao.PrintSpec();

	//prtoected.Calc(paras_file_volup);
	
	//cout << "\nprtoected.Calc(paras_file);" << endl;
	//prtoected.PrintResult();

	///ao.Simulation2(paras_file, 100, true);
	vector<double> apath;
	apath = get_a_path_from_csv("kith.csv");
	//notprot20180828.Simulation3(paras_20180828, apath, true);
	//notprot20180828.Simulation2(paras_20180828, 5000, true);
	notprot.CalcMC(paras_file, 10000);
	notprot.PrintResult();
	notprot.CalcMC_discrete(paras_file_div, 10000);
	notprot.PrintResult();

	//ao.Simulation3(paras_file_volup, apath, true);

	//ao.Simulation2(paras_file, 10, true);
}
void test_autocall_final()
{
	MarketParameters paras;
	MarketParameters paras_volup;
	paras = set_paras();
	paras_volup = set_paras();
	
	double refprice = 297.22;
	double putstrike = refprice;
	double putstrike90 = refprice*0.9;

	int hitflag = 0;
	signed int exd_3y = 44436; //3y 
	int nb_autocall = 6;
	signed int auto_date[7] = { -1,43524,43705,43889,44071,44255,44436 };
	double auto_strike[7] = { -1,297.22,297.22,282.359,282.359,267.498,267.498 };
	double auto_strike_hifive[7] = { -1,297.22,297.22,297.22,297.22,297.22,297.22 };
	double auto_strike959085[7] = { -1,282.359,282.359,267.498,267.498,252.637,252.637 };
	double auto_strike_shift[7] = { -1,297.22,297.22,230.3455,230.3455,230.3455,230.3455 };
	double auto_coupon[7] = { -1, 0.0230,	0.0460,	0.0690,	0.0920,	0.1150,	0.1380 };
	double auto_coupon32[7] = { -1, 0.0160,	0.0320,	0.0480,	0.0640,	0.0800,	0.0960 };
	double auto_coupon28[7] = { -1, 0.0140,	0.0280,	0.0420,	0.0560,	0.0700,	0.0840 };

	double auto_ki_barrier = refprice*0.6;
	double auto_ki_barrier70 = refprice * 0.7;
	double auto_ki_barrier80 = refprice * 0.8;

	double auto_ki_barrier_shift = refprice*0.6;
	double auto_dummy_coupon = auto_coupon[6];
	double auto_dummy_coupon32 = auto_coupon32[6];

	double auto_put_strike = refprice*1.0; // if put_strike=0, notional protected
	double put_strike_notional_protect = refprice*0.0;

	PayoffAutocallStd autoPayoff(nb_autocall, auto_date, auto_strike, auto_coupon, auto_ki_barrier, auto_put_strike, auto_dummy_coupon, refprice);
	PayoffAutocallStd autoPayoff_shift(nb_autocall, auto_date, auto_strike_shift, auto_coupon, auto_ki_barrier, auto_put_strike, auto_dummy_coupon, refprice);
	PayoffAutocallStd autoPayoff_KI_shift(nb_autocall, auto_date, auto_strike, auto_coupon, auto_ki_barrier_shift, auto_put_strike, auto_dummy_coupon, refprice);
	PayoffAutocallStd autoPayoff959085(nb_autocall, auto_date, auto_strike959085, auto_coupon, auto_ki_barrier, auto_put_strike, auto_dummy_coupon, refprice);
	PayoffAutocallStd autoPayoff_hifive(nb_autocall, auto_date, auto_strike_hifive, auto_coupon, auto_ki_barrier70, put_strike_notional_protect, auto_dummy_coupon, refprice);
	PayoffAutocallStd autoPayoff_notional_protect(nb_autocall, auto_date, auto_strike_hifive, auto_coupon28, auto_ki_barrier80, put_strike_notional_protect, auto_coupon28[6], refprice);

	AutocallOption AutoKOSPI(refprice, exd_3y, autoPayoff, hitflag);
	AutocallOption AutoKOSPIMC(refprice, exd_3y, autoPayoff, hitflag);

	AutocallOption AutoKOSPI_hit(refprice, exd_3y, autoPayoff, 1);
	AutocallOption AutoKOSPI959085(refprice, exd_3y, autoPayoff959085, hitflag);
	AutocallOption AutoKOSPI_hifive(refprice, exd_3y, autoPayoff_hifive, 0);
	AutocallOption AutoKOSPI_hifive_hit(refprice, exd_3y, autoPayoff_hifive, 1);
	AutocallOption AutoKOSPI_shift(refprice, exd_3y, autoPayoff_shift, hitflag);
	AutocallOption AutoKOSPI_KI_shift(refprice, exd_3y, autoPayoff_KI_shift, hitflag);
	AutocallOption AutoKOSPI_notional_protect(refprice, exd_3y, autoPayoff_notional_protect, hitflag);
	AutocallOption AutoKOSPI_notional_protect_flatvol(refprice, exd_3y, autoPayoff_notional_protect, hitflag);
	AutocallOption AutoKOSPI_mc(refprice, exd_3y, autoPayoff, hitflag);
	AutocallOption AutoKOSPI_KI(refprice, exd_3y, autoPayoff, 1);
	AutocallOption AutoKOSPI_mc_KI(refprice, exd_3y, autoPayoff, 1);

	AutoKOSPI.Calc(paras);
	cout << "\nAutoKOSPI.Calc(paras);" << endl;
	AutoKOSPI.PrintResult();

	AutoKOSPI.CalcMC(paras,10000);
	cout << "\nAutoKOSPI.CalcMC(paras);" << endl;
	AutoKOSPI.PrintResult();

}



MarketParameters init_paras_file(const char* vol_csv, const char* rate_csv, const char* div_csv)
{
	MarketParameters paras;

	vector<double> vol_term;
	vector<double> vol_strike;
	vector<vector<double> > Ivol;
	vector<double> rate;
	vector<double> rate_term;
	vector<double> qrate;
	vector<double> qrate_term;


	ifstream infile(vol_csv); // for example
	string line = "";
	string word = "";
	vector<vector<double> >rows;
	getline(infile, line); //first vol line, strikes
	stringstream strstr(line);
	
	for(int i=0;i<5;i++) //skip first 5 cells
		getline(strstr,word,',');
	while (getline(strstr, word, ',')) {
		vol_strike.push_back(stod(word));
	}
	
	getline(infile, line); //skip the second line

	while (getline(infile, line)) {
		stringstream strstr(line);
		string word = "";
		
		//cells(3,1), cells(4,1).... => vol_term
		getline(strstr, word, ',');
		vol_term.push_back(stod(word));

		for (int i = 0; i<4; i++) //skip first 4 cells
			getline(strstr, word, ',');

		vector<double> avol;
		while (getline(strstr, word, ',')) {
			avol.push_back(stod(word)/100.0);
		}
		Ivol.push_back(avol);
	}

	Vol vol(vol_term, vol_strike, Ivol);
	infile.close();

	ifstream infile_rate(rate_csv);
	//line = ""; //not needed?
	getline(infile_rate, line); //first line, rate term
	//stringstream().swap(strstr);
	strstr = stringstream(line);

	//word = "";
	getline(strstr, word, ','); //skip cells(1,1)
	while (getline(strstr, word, ',')) {
		rate_term.push_back(stod(word));
	}
	getline(infile_rate, line); //skip 2nd line
	getline(infile_rate, line); //3rd line
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip cells(3,1)
	while (getline(strstr, word, ',')) {
		rate.push_back(stod(word)/100.0);
	}
	Rate r(rate, rate_term);
	infile_rate.close();

	ifstream infile_qrate(div_csv);
	//line = ""; //not needed?
	getline(infile_qrate, line); //skip first line
	getline(infile_qrate, line); //skip 2nd line
	getline(infile_qrate, line); //3rd line qrate_Term
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip cells(3,1)
	while (getline(strstr, word, ',')) {
		qrate_term.push_back(stod(word));
	}

	getline(infile_qrate, line); //line4
	strstr = stringstream(line);
	getline(strstr, word, ','); 
	double spot = stod(word);
	while (getline(strstr, word, ',')) {
		//skip 
	}

	getline(infile_qrate, line); //line5
	strstr = stringstream(line);
	getline(strstr, word, ',');
	signed int vd = stoul(word);
	while (getline(strstr, word, ',')) {
		qrate.push_back(stod(word)/100.0);
	}
	Rate q(qrate, qrate_term);
	infile_qrate.close();

	return MarketParameters(vd, spot, vol, r, q);
	//return MarketParameters(vol,r,q);
}

MarketParameters init_paras_file(const char* vol_csv, const char* rate_csv, const char* div_csv, const char* div_discrete_csv)
{
	MarketParameters paras;

	vector<double> vol_term;
	vector<double> vol_strike;
	vector<vector<double> > Ivol;
	vector<double> rate;
	vector<double> rate_term;
	vector<double> qrate;
	vector<double> qrate_term;
	vector<signed int> exdate;
	vector<double> amount;

	ifstream infile(vol_csv); // for example
	string line = "";
	string word = "";
	vector<vector<double> >rows;
	getline(infile, line); //first vol line, strikes
	stringstream strstr(line);

	for (int i = 0; i<5; i++) //skip first 5 cells
		getline(strstr, word, ',');
	while (getline(strstr, word, ',')) {
		vol_strike.push_back(stod(word));
	}

	getline(infile, line); //skip the second line

	while (getline(infile, line)) {
		stringstream strstr(line);
		string word = "";

		//cells(3,1), cells(4,1).... => vol_term
		getline(strstr, word, ',');
		vol_term.push_back(stod(word));

		for (int i = 0; i<4; i++) //skip first 4 cells
			getline(strstr, word, ',');

		vector<double> avol;
		while (getline(strstr, word, ',')) {
			avol.push_back(stod(word) / 100.0);
		}
		Ivol.push_back(avol);
	}

	Vol vol(vol_term, vol_strike, Ivol);
	infile.close();

	ifstream infile_rate(rate_csv);
	//line = ""; //not needed?
	getline(infile_rate, line); //first line, rate term
								//stringstream().swap(strstr);
	strstr = stringstream(line);

	//word = "";
	getline(strstr, word, ','); //skip cells(1,1)
	while (getline(strstr, word, ',')) {
		rate_term.push_back(stod(word));
	}
	getline(infile_rate, line); //skip 2nd line
	getline(infile_rate, line); //3rd line
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip cells(3,1)
	while (getline(strstr, word, ',')) {
		rate.push_back(stod(word) / 100.0);
	}
	Rate r(rate, rate_term);
	infile_rate.close();

	//continuous dividend
	ifstream infile_qrate(div_csv);
	getline(infile_qrate, line); //skip first line
	getline(infile_qrate, line); //skip 2nd line
	getline(infile_qrate, line); //3rd line qrate_Term
	strstr = stringstream(line);
	getline(strstr, word, ','); //skip cells(3,1)
	while (getline(strstr, word, ',')) {
		qrate_term.push_back(stod(word));
	}

	getline(infile_qrate, line); //line4
	strstr = stringstream(line);
	getline(strstr, word, ',');
	double spot = stod(word);
	while (getline(strstr, word, ',')) {
		//skip 
	}

	getline(infile_qrate, line); //line5
	strstr = stringstream(line);
	getline(strstr, word, ',');
	signed int vd = stoul(word);
	while (getline(strstr, word, ',')) {
		qrate.push_back(stod(word) / 100.0);
	}
	Rate q(qrate, qrate_term);
	infile_qrate.close();
	
	//discrete Dividend
	infile= ifstream(div_discrete_csv);
	while (getline(infile, line)) {
		strstr = stringstream(line);
		getline(strstr, word, ','); //col=1
		exdate.push_back(stoul(word));
		getline(strstr, word, ','); //col=2
		amount.push_back(stod(word));
	}
	Dividend div(amount,exdate);
	infile.close();

	return MarketParameters(vd, spot, vol,r,q,div);
	
}
int main()
{
	using clock = std::chrono::system_clock;
	using sec = std::chrono::duration<double>;
	auto before = clock::now();
	
	//sec duration = clock::now() - before;
	//cout << "\n"<<duration.count() << "s" << std::endl;
	//MarketParameters p = set_file("volcsv.csv");;
	
	//test_vanilla_final();
	//test_vanilla_iofile();
	
	//test_autocall_final();
	test_autocall_final_iofile();
}

