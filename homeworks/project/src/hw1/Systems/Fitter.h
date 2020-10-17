#pragma once
#include <functional>
#include <vector>
#include <UGM/point.h>

class Fitter
{
public:

	void set_data(std::vector<Ubpa::pointf2 >& x)
	{
		data = x;
		fit_data();
	}

	double f(double x) const { return  rst(x); }

	virtual void fit_data() = 0;

protected:
	std::vector<Ubpa::pointf2> data;
	std::function<double(double)> rst;

	double x(int index) { return data[index][0]; }
	double y(int index) { return data[index][1]; }
};

class PolynomialFitter :public Fitter
{
public:
	void fit_data() override;
private:
	double diff(double y2, double y1, double x2, double x1) { return (y2 - y1) / (x2 - x1); }
};

class GaussianFitter :public Fitter
{
public:
	void fit_data() override;
};

class LSFitter :public Fitter
{
public:
	void fit_data() override;
};

class MLSFitter : public Fitter
{
public:
	void fit_data() override;
};

class RidgeFitter :public Fitter
{
public:
	void fit_data() override;
};
