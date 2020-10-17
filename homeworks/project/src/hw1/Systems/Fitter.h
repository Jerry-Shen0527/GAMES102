#pragma once
#include <functional>
#include <vector>
#include <UGM/point.h>
using Ubpa::pointf2;
class Fitter
{
public:

	void set_data(std::vector<pointf2 >& x)//fit_data will be automatically called.
	{
		data = x;
		fit_data();
	}
	//comfortably extracting f(x) from the fitted model.
	double f(double x) const
	{
		if (rst == nullptr)
		{
			return 0;
		}
		return  rst(x);
	}

	bool ready() { return rst != nullptr; }

protected:
	std::vector<pointf2> data;
	std::function<double(double)> rst;//Use a closure to store the result
	virtual void fit_data() = 0;

	double x(int index) { return data[index][0]; }
	double y(int index) { return data[index][1]; }

	bool fitted = false;
};

class PolynomialFitter :public Fitter
{
public:
	PolynomialFitter() {}
	void fit_data() override;
private:
	static double diff(double y2, double y1, double x2, double x1) { return (y2 - y1) / (x2 - x1); }
};

class GaussianFitter :public Fitter
{

public:
	GaussianFitter(){}
	GaussianFitter(double sigma):sigma(sigma){}
	void fit_data() override;
private:
	double gauss(double x, double x0, double sigma)
	{
		return exp(-(x - x0) * (x - x0) / sigma / sigma);
	}
	double sigma=100;
};

class LSFitter :public Fitter
{
public:
	LSFitter() {}
	LSFitter(int n) :max_n(n) {}
	void fit_data() override;

private:
	int max_n = 4;
};

class RidgeFitter :public Fitter
{
public:
	RidgeFitter() {}
	RidgeFitter(int n, double alpha_) :max_n(n), alpha(alpha_) {}
	void fit_data() override;
private:
	int max_n = 4;
	double alpha = 0.2;
};
