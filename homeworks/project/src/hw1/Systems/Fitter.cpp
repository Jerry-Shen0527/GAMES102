#include "Fitter.h"

#include <iostream>

#include "eigen3/Eigen/Eigen"

void PolynomialFitter::fit_data()
{
	using namespace std;
	if (data.empty())
	{
		return;
	}

	vector<double> diff_quotion;
	vector<double> new_diff_quotion;
	vector<double> rst_diff_quotion;
	diff_quotion.resize(data.size());
	new_diff_quotion.resize(data.size());
	rst_diff_quotion.resize(data.size());

	for (int i = 0; i < data.size(); ++i)
	{
		new_diff_quotion[0] = y(i);
		for (int j = 1; j <= i; ++j)
		{
			new_diff_quotion[j] = diff(new_diff_quotion[j - 1], diff_quotion[j - 1], x(i), x(i - j));
		}
		std::swap(diff_quotion, new_diff_quotion);
		rst_diff_quotion[i] = diff_quotion[i];
	}

	rst = [rst_diff_quotion, this](const double x)->double
	{
		double rst = 0;
		double temp = 1;

		for (int i = 0; i < data.size() - 1; ++i)
		{
			rst += rst_diff_quotion[i] * temp;
			temp *= x - this->x(i);
		}

		return rst + rst_diff_quotion[data.size() - 1] * temp;
	};
}

void GaussianFitter::fit_data()
{
	using namespace Eigen;

	auto size = data.size();

	SparseMatrix<double> mat(size, size);
	VectorXd right_hand(size);

	std::vector<Triplet<double>> triplets;

	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			triplets.emplace_back(i, j, gauss(x(i), x(j), sigma));
		}
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());

	for (int i = 0; i < size; ++i)
	{
		right_hand(i) = y(i);
	}

	SparseLU<SparseMatrix<double>> sparse_lu;
	sparse_lu.compute(mat);
	if (sparse_lu.info() != Eigen::Success)
	{
		std::cerr << "unsolved" << std::endl;
		return;
	}
	VectorXd solution = sparse_lu.solve(right_hand);

	rst = [solution, this](double x)->double
	{
		double result = 0;
		for (int i = 0; i < data.size(); ++i)
		{
			result += solution(i) * gauss(x, this->x(i), sigma);
		}
		return result;
	};
}

void LSFitter::fit_data()
{
	using namespace std;

	auto size = data.size();
	vector<double> coeffi;
	coeffi.resize(data.size());

	using namespace Eigen;
	VectorXd target(size);

	for (int i = 0; i < size; ++i)
	{
		target(i) = y(i);
	}

	SparseMatrix<double> mat(size, max_n + 1);

	vector<Triplet<double>> triplets;
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j <= max_n; ++j)
		{
			triplets.emplace_back(i, j, pow(x(i), j));
		}
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());
	VectorXd right_hand_side = mat.transpose() * target;
	auto mat_to_be_solved = mat.transpose() * mat;
	SparseLU<SparseMatrix<double>> sparse_lu;
	sparse_lu.compute(mat_to_be_solved);
	if (sparse_lu.info() != Eigen::Success)
	{
		std::cerr << "unsolved" << std::endl;
		return;
	}
	VectorXd solution = sparse_lu.solve(right_hand_side);

	rst = [solution, this](double x)->double
	{
		double result = 0;
		double temp = 1;
		for (int i = 0; i <= max_n; ++i)
		{
			result += temp * solution(i);
			temp *= x;
		}
		return result;
	};
}

void RidgeFitter::fit_data()
{
	using namespace std;

	auto size = data.size();
	vector<double> coeffi;
	coeffi.resize(data.size());

	using namespace Eigen;
	VectorXd target(size);

	for (int i = 0; i < size; ++i)
	{
		target(i) = y(i);
	}

	SparseMatrix<double> mat(size, max_n + 1);

	vector<Triplet<double>> triplets;
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j <= max_n; ++j)
		{
			triplets.emplace_back(i, j, pow(x(i), j));
		}
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());
	VectorXd right_hand_side = mat.transpose() * target;
	auto mat_to_be_solved = mat.transpose() * mat + alpha * MatrixXd::Identity(max_n + 1, max_n + 1);
	SparseLU<SparseMatrix<double>> sparse_lu;
	sparse_lu.compute(mat_to_be_solved);
	if (sparse_lu.info() != Eigen::Success)
	{
		std::cerr << "unsolved" << std::endl;
		return;
	}
	VectorXd solution = sparse_lu.solve(right_hand_side);

	rst = [solution, this](double x)->double
	{
		cout << solution;
		double result = 0;
		double temp = 1;
		for (int i = 0; i <= max_n; ++i)
		{
			result += temp * solution(i);
			temp *= x;
		}
		return result;
	};
}