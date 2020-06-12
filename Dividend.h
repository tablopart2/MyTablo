#pragma once
#include <vector>
#include "Rate.h"

using namespace std;

class Dividend {
public:
	Dividend() {}
	Dividend(vector<double> _amount, vector<signed int> _exdate) : amount(_amount), exdate(_exdate){};
	double getLumpsum(signed int start, signed int end) const;
	Rate getDividendCurve(std::vector<double>& rts, signed int vd, double spot);
private:
	vector<double> amount;
	vector<signed int> exdate;
};