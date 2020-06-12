#include "Dividend.h"
#include <algorithm> //std::find
#include <vector>
#include <iostream>

using namespace std;

double Dividend::getLumpsum(signed int start_date, signed int end_date) const
{
	double lumpsum = 0;

	for (signed int t = start_date; t <= end_date; t++) {
		vector<signed int>::const_iterator found = find(exdate.begin(), exdate.end(), t);
		if (found != exdate.end()) {
			auto index = distance(exdate.begin(), found);
			lumpsum += amount[index];
		}
	}

	return lumpsum;
}

Rate Dividend::getDividendCurve(vector<double>&rts, signed int vd, double spot)
{
	vector<double> r;
	//vector<double>::iterator it = rts.begin();
	for (vector<double>::iterator it = rts.begin(); it != rts.end(); it++) {
		double rat = getLumpsum(vd + 1, vd + int((*it) * 365)) / spot / (*it);
		//cout << rat << endl;
		r.push_back(rat);
	}
	
	return Rate(r, rts);
}
