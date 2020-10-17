#include "Fitter.h"

void PolynomialFitter::fit_data()
{
	using namespace std;
	if (data.empty())
	{
		return;
	}

	vector<double> diff_quotion;

	double temp = 0;
	for (int i = 0; i < data.size(); ++i)
	{
		temp = y(i);
		for (int j = 0; j < i; ++j)
		{
			diff_quotion[j] = temp;
			temp = diff_quotion[j];

			diff_quotion[j] = temp;
		}
		diff_quotion.push_back(temp);
	}

	rst =
		[&diff_quotion, this](double)->double
	{
		return  0;
	};
}