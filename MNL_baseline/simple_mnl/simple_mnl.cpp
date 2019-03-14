// Actually a logit exercise
// Works with example solver
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#define n_times 4
#define n_options 7
#define n_samples 14

namespace
{
using CppAD::AD;

int data[n_samples][4];

template <typename T>
void printMatrix(const char *text, T mat, std::size_t N, std::size_t M, int width)
{
	std::cout << text << std::endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
			std::cout << std::setw(width) << std::fixed << std::setprecision(2)
					  << mat[i][j] << std::setw(width);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// helper function for printing vectors (debugging)
template <typename T>
void printVector(const char *text, T mat, std::size_t N, int width)
{
	std::cout << text << std::endl;

	for (int i = 0; i < N; i++)
	{
		std::cout << std::setw(width) << std::fixed << std::setprecision(2)
				  << mat[i] << std::setw(width);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


class FG_eval
{
  public:
	typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
	void operator()(ADvector &fg, const ADvector &x)
	{
		assert(fg.size() == 1);
		assert(x.size() == 8);

		//data (normally would be imported, but this time hand generated)
		//populate cost and price matrices
		//data needs to be sorted by time, then by option ID
		int prices[n_times][n_options] = {{0}}; //price matrix
		int counts[n_times][n_options] = {{0}}; //count matrix
		int t;
		int i;
		printMatrix("DATA: ", data, n_samples, 4, 6);
		for (int s = 0; s < n_samples; s++)
		{
			t = data[s][0];
			i = data[s][1];
			prices[t][i] = data[s][2];
			counts[t][i] = data[s][3];
		}
		printMatrix("PRICES: ", prices, n_times, n_options, 3);
		printMatrix("COUNTS: ", counts, n_times, n_options, 3);
		printVector("X:", x, 7, 3);

		//calculating utility matrix
		AD<double> utils[n_times][n_options] = {{0}};
		for (int t = 0; t < n_times; t++)
		{
			for (int i = 0; i < n_options; i++)
			{
				//std::cout << x[i] << std::endl;
				if (counts[t][i] > 0)
				{
					utils[t][i] = x[i] + x[7] * prices[t][i];
				}
			}
		}
		printMatrix("UTILS: ", utils, n_times, n_options, 3);

		//calculating prob matrix
		//first need to calculate avg prob per time period
		AD<double> avg_p[n_times]; // array avg prob over all options in a time period
		for (int t = 0; t < n_times; t++)
		{
			AD<double> sumprob = 0;
			for (int i = 0; i < n_options; i++)
			{
				sumprob += exp(utils[t][i]);
			}
			avg_p[t] = sumprob;
		}
		
		//calculating matrix of probabilities
		AD<double> probs[n_times][n_options] = {{0}};
		for (int t = 0; t < n_times; t++)
		{
			for (int i = 0; i < n_options; i++)
			{
				probs[t][i] = exp(utils[t][i]) / avg_p[t];
			}
		}

		// printMatrix("PROBS: ", probs, n_times, n_options, 3);

		//building LL
		for (int t = 0; t < n_times; t++)
		{
			for (int i = 0; i < n_options; i++)
			{
				// f(x)
				fg[0] -= counts[t][i] * log(probs[t][i]);
			}
		}

		return;
	}
};
} // namespace

int optimize()
{
	bool ok = true;

	size_t i;
	typedef CPPAD_TESTVECTOR(double) Dvector;

	// number of independent variables (domain dimension for f and g)
	size_t nx = 8;
	// number of constraints (range dimension for g)
	size_t ng = 0;
	// initial value of the independent variables
	Dvector xi(nx);
	xi[0] = 1.0;
	xi[1] = 1.0;
	xi[2] = 1.0;
	xi[3] = 1.0;
	xi[4] = 1.0;
	xi[5] = 1.0;
	xi[6] = 1.0;
	xi[7] = 1.0;
	// lower and upper limits for x
	Dvector xl(nx), xu(nx);
	for (i = 0; i < nx; i++)
	{
		xl[i] = -1e19;
		xu[i] = 1e19;
	}
	// lower and upper limits for g
	Dvector gl(ng), gu(ng);

	// object that computes objective and constraints
	FG_eval fg_eval;

	// options
	std::string options;
	// turn off any printing
	//options += "Integer print_level  0\n";
	//options += "String  sb           yes\n";
	// maximum number of iterations
	options += "Integer max_iter     200\n";
	// approximate accuracy in first order necessary conditions;
	// see Mathematical Programming, Volume 106, Number 1,
	// Pages 25-57, Equation (6)
	options += "Numeric tol          1e-8\n";
	// derivative testing
	options += "String  derivative_test            second-order\n";
	// maximum amount of random pertubation; e.g.,
	// when evaluation finite diff
	options += "Numeric point_perturbation_radius  4.\n";

	// place to return solution
	CppAD::ipopt::solve_result<Dvector> solution;

	// solve the problem
	CppAD::ipopt::solve<Dvector, FG_eval>(
		options, xi, xl, xu, gl, gu, fg_eval, solution);
	//
	// Check some of the solution values
	//
	ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
	//

	return ok;
}

int main()
{
	int ok;

	// initial data
	int data_filler[n_samples][4] =
		{
			{0, 1, 2, 10},
			{0, 2, 4, 2},
			{0, 3, 6, 4},
			{0, 5, 10, 11},
			{0, 6, 12, 23},
			{1, 2, 4, 8},
			{1, 3, 6, 3},
			{1, 4, 8, 2},
			{1, 6, 12, 1},
			{2, 1, 2, 6},
			{2, 2, 4, 32},
			{3, 1, 2, 66},
			{3, 3, 6, 21},
			{3, 5, 10, 19}};

	for (int s = 0; s < n_samples; s++)
	{
		for (int i = 0; i < 4; i++)
		{
			data[s][i] = data_filler[s][i];
		}
	}

	// test multipler solves with different data
	for (int j = 0; j < 5; j++)
	{
		ok = optimize();
		for (int s = 0; s < n_samples; s++)
		{
			data[s][3] = rand() % 100;
		}
	}

	return ok;
}
