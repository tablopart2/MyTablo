#ifndef RATE_H
#define RATE_H

#include <vector>

class Rate{
public:
	Rate(){}
	Rate(std::vector<double> _r, std::vector<double> _ts):r(_r),ts(_ts){};
	void set_const_rate(double v); //initialize value with constant
	double getForward(double t) const;
	double getIntpRate(double t) const;
	std::vector<double> get_rate() const;
	void print() const;
	int getSize() const;

private:
	std::vector<double> r;
	std::vector<double> ts;

};
#endif