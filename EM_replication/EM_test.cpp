// Actually a logit exercise
// Works with example solver
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <math.h>
#include "csv.h"

#define n_times 245 // number of
#define n_options 8 //
#define n_types 14  // number of customer types

namespace
{
using CppAD::AD;

// helper function for printing matrices (debugging)
template <typename T>
void printMatrix(T mat, std::size_t N, std::size_t M, int width)
{
	std::cout << "\n Printing Matrix : \n";

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
			std::cout << std::setw(width) << mat[i][j] << std::setw(width);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// helper function for printing vectors (debugging)
template <typename T>
void printVector(T mat, std::size_t N, int width)
{
	std::cout << "\n Printing Vector : \n";

	for (int i = 0; i < N; i++)
	{
		std::cout << std::setw(width) << mat[i] << std::setw(width);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// import preference matrix
// returns pref_matrix: a n_types x (n_options + 1) matrix
AD<int> **import_prefs(const char *pref_filename)
{
	io::CSVReader<10> in(pref_filename);
	in.read_header(io::ignore_no_column, "cust_type", "rank_1", "rank_2", "rank_3",
				   "rank_4", "rank_5", "rank_6", "rank_7", "rank_8", "rank_9");
	int cust_type;
	int pref_rank[n_options + 1];
	int row_counter = 0;
	AD<int> **pref_matrix = 0;
	pref_matrix = new AD<int> *[n_types];
	while (in.read_row(cust_type, pref_rank[0], pref_rank[1], pref_rank[2],
					   pref_rank[3], pref_rank[4], pref_rank[5],
					   pref_rank[6], pref_rank[7], pref_rank[8]))
	{
		pref_matrix[row_counter] = new AD<int>[n_options + 1];
		for (int i = 0; i < n_options + 1; i++)
		{
			pref_matrix[row_counter][i] = pref_rank[i];
		}
		row_counter++;
	}
	return pref_matrix;
}

// import availability matrix
// returns avail_matrix: a n_times x n_options matrix
AD<int> **import_availability(const char *avail_filename)
{
	io::CSVReader<9> in(avail_filename);
	in.read_header(io::ignore_no_column, "T", "prod_1", "prod_2", "prod_3",
				   "prod_4", "prod_5", "prod_6", "prod_7", "prod_8");
	int t;
	int prod[n_options];
	int row_counter = 0;
	AD<int> **avail_matrix = 0;
	avail_matrix = new AD<int> *[n_times];
	while (in.read_row(t, prod[0], prod[1], prod[2], prod[3], prod[4],
					   prod[5], prod[6], prod[7]))
	{
		avail_matrix[row_counter] = new AD<int>[n_options];
		for (int i = 0; i < n_options; i++)
		{
			avail_matrix[row_counter][i] = prod[i];
		}
		row_counter++;
	}
	return avail_matrix;
}

// import transaction vector
// returns trans_vector: a length T vector
AD<int> *import_transactions(const char *trans_filename)
{
	io::CSVReader<2> in(trans_filename);
	in.read_header(io::ignore_no_column, "T", "prod_num");
	int t;
	int prod_num;
	int row_counter = 0;
	AD<int> *trans_vec = new AD<int>[n_times];
	while (in.read_row(t, prod_num))
	{
		trans_vec[row_counter] = prod_num;
		row_counter++;
	}
	return trans_vec;
}

// build compatible type (mu_t) sets
// returns mu_matrix: a n_times x n_types matrix
// each col is a type, each row is a time period
// if value is 1, then that type is compatible in that time period
AD<int> **build_mu_mat(AD<int> pref_matrix[n_types][n_options + 1],
					   AD<int> avail_matrix[n_times][n_options], AD<int> trans_vec[n_times])
{
	printMatrix(pref_matrix, 10, n_options + 1, 3);
	printMatrix(avail_matrix, 10, n_options, 3);
	printVector(trans_vec, 10, 3);

	AD<int> **mu_matrix = 0;
	mu_matrix = new AD<int> *[n_times];
	for (int t = 0; t < n_times; t++)
	{
		mu_matrix[t] = new AD<int>[n_types];
		for (int i = 0; i < n_types; i++)
		{
			avail_matrix[t][i] = 0;
		}
	}
	return mu_matrix;
}

// Problem formulated here
class FG_eval
{
  public:
	typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
	void operator()(ADvector &fg, const ADvector &x)
	{
		//Compatability code for ipopt ##################################################
		assert(fg.size() == 1);
		assert(x.size() == 1);
		AD<double> x1 = x[0];
		fg[0] = 1;
		//Compatability code for ipopt ##################################################

		AD<int> **pref_matrix = import_prefs("data/hotel_5/PrefListsBuyUpH5.csv");
		AD<int> **avail_matrix = import_availability("data/hotel_5/AvailabilityH5.csv");
		AD<int> *trans_vec = import_transactions("data/hotel_5/TransactionsH5.csv");

		AD<int> **mu_matrix = build_mu_mat(pref_matrix, avail_matrix, trans_vec);

		// Debugging prints
		std::cout << "PREFERENCE MATRIX" << std::endl;
		printMatrix(pref_matrix, 10, n_options + 1, 3);
		std::cout << "AVAILABILITY MATRIX" << std::endl;
		printMatrix(avail_matrix, 10, n_options, 3);
		std::cout << "TRANSACTION MATRIX" << std::endl;
		printVector(trans_vec, 10, 3);

		return;
	}
};
} // namespace

int main()
{
	//Compatability code for ipopt ##################################################
	bool ok = true;
	size_t i;
	typedef CPPAD_TESTVECTOR(double) Dvector;

	// number of independent variables (domain dimension for f and g)
	size_t nx = 1;
	// number of constraints (range dimension for g)
	size_t ng = 0;
	// initial value of the independent variables
	Dvector xi(nx);
	xi[0] = 1.0;

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
	//Compatability code for ipopt ##################################################

	std::cout << "TEST DONE" << std::endl;
}