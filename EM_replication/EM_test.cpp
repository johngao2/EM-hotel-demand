// Actually a logit exercise
// Works with example solver
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <cppad/ipopt/solve.hpp>
#include "csv.h"

#define n_times 245 // number of time steps
#define n_options 8 // number of products
#define n_types 14  // number of customer types

namespace
{
using CppAD::AD;

// helper function for printing matrices (debugging)
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

// import preference matrix (sigma_matrix)
// this is equivalent to the sigma table where sigma_i(j_t)
// = rank of j_t in customer type i
// each row is a cust type, each col is a product, each value is a rank
// first column represents non-purchase option
// 0 is most preferred, n_options+1 is least preferred
// returns sigma_matrix: a n_types x (n_options + 1) matrix
int **import_prefs(const char *pref_filename)
{
	// import original preference matrix from data
	io::CSVReader<10> in(pref_filename);
	in.read_header(io::ignore_no_column, "cust_type", "rank_1", "rank_2", "rank_3",
				   "rank_4", "rank_5", "rank_6", "rank_7", "rank_8", "rank_9");
	int cust_type;
	int pref_rank[n_options + 1];
	int row_counter = 0;
	int **pref_matrix = 0;
	pref_matrix = new int *[n_types];
	while (in.read_row(cust_type, pref_rank[0], pref_rank[1], pref_rank[2],
					   pref_rank[3], pref_rank[4], pref_rank[5],
					   pref_rank[6], pref_rank[7], pref_rank[8]))
	{
		pref_matrix[row_counter] = new int[n_options + 1];
		for (int i = 0; i < n_options + 1; i++)
		{
			pref_matrix[row_counter][i] = pref_rank[i];
		}
		row_counter++;
	}

	// convert to sigma matrix
	int **sigma_matrix = 0;
	sigma_matrix = new int *[n_types];
	bool viable; // tracks if an option is viable (i.e. ranked lower than nonpurchase)
	// constructs sigma matrix for viable options
	for (int i = 0; i < n_types; i++)
	{
		sigma_matrix[i] = new int[n_options + 1];
		viable = 1;
		for (int rank = 0; rank < n_options + 1; rank++)
		{
			int k = pref_matrix[i][rank];
			if (viable == 1)
			{
				if (k == 0)
				{
					viable = 0;
					sigma_matrix[i][k] = rank + 1;
				}
				else
				{
					sigma_matrix[i][k] = rank + 1;
				}
			}
		}
	}

	// zero index rankings, and fill nonviable options with rank (n_options+1)
	for (int i = 0; i < n_types; i++)
	{
		for (int k = 0; k < n_options + 1; k++)
		{
			if (sigma_matrix[i][k] == 0)
			{
				sigma_matrix[i][k] = n_options + 1;
			}
			else
			{
				sigma_matrix[i][k]--;
			}
		}
	}

	return sigma_matrix;
}

// import availability matrix
// returns avail_matrix: a n_times x n_options matrix
int **import_availability(const char *avail_filename)
{
	io::CSVReader<9> in(avail_filename);
	in.read_header(io::ignore_no_column, "T", "prod_1", "prod_2", "prod_3",
				   "prod_4", "prod_5", "prod_6", "prod_7", "prod_8");
	int t;
	int prod[n_options];
	int row_counter = 0;
	int **avail_matrix = 0;
	avail_matrix = new int *[n_times];
	while (in.read_row(t, prod[0], prod[1], prod[2], prod[3], prod[4],
					   prod[5], prod[6], prod[7]))
	{
		avail_matrix[row_counter] = new int[n_options];
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
int *import_transactions(const char *trans_filename)
{
	io::CSVReader<2> in(trans_filename);
	in.read_header(io::ignore_no_column, "T", "prod_num");
	int t;
	int prod_num;
	int row_counter = 0;
	int *trans_vec = new int[n_times];
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
int **build_mu_mat(int **sigma_matrix, int **avail_matrix, int *trans_vec)
{
	int **mu_matrix = 0;
	mu_matrix = new int *[n_times];
	int *ranking;
	bool compatible; // boolean var to keep track of if a type is compatible

	for (int t = 0; t < n_times; t++)
	{
		mu_matrix[t] = new int[n_types];
		int j_t = trans_vec[t];			  // subtracted for zero indexing
		for (int i = 0; i < n_types; i++) // iterate over each cust type
		{
			compatible = 1;
			ranking = sigma_matrix[i];
			for (int k = 1; k < n_options + 1; k++) // iterate over each option
			{
				bool available = avail_matrix[t][k];
				if (available == 1 && k != j_t) // if item is available and not selected
				{
					if (sigma_matrix[i][k] < sigma_matrix[i][j_t])
					// if any item is ranked better than chosen item
					{
						compatible = 0;
					}
				}
			}
			// check if nonpurchase is preferred to chosen item
			if (sigma_matrix[i][0] < sigma_matrix[i][j_t])
			{
				compatible = 0;
			}
			mu_matrix[t][i] = compatible;
		}
	}
	return mu_matrix;
}

// initialize m vector
// returns m_vec: length N
double *initialize_m_vec()
{
	double *m_vec = new double[n_times];
	for (int i = 0; i < n_types; i++)
	{
		m_vec[i] = 0;
	}
	return m_vec;
}

// build p_sigma matrix (calculates customer type probabilities
// based on the data and previous guess for xi)
// p_sigma is confusingly denoted as x_it in the pseudocode
double **build_cust_type_probs(double *x, int **mu_matrix)
{
	double **p_sigma_matrix = 0;
	p_sigma_matrix = new double *[n_times];
	for (int t = 0; t < n_times; t++)
	{
		p_sigma_matrix[t] = new double[n_types];

		// first calculate sum probabilities in each time
		// this is the denominator sum(x_h) for all compatible types h
		double sumprob = 0;
		for (int i = 0; i < n_types; i++)
		{
			sumprob += mu_matrix[t][i] * x[i];
		}

		// update type probability if that type is compatible
		for (int i = 0; i < n_types; i++)
		{
			if (mu_matrix[t][i] == 1) // set was compatible, so update this value
			{
				p_sigma_matrix[t][i] = x[i] / sumprob;
			}
			else
			{
				p_sigma_matrix[t][i] = 0;
			}
		}
	}
	return p_sigma_matrix;
}

// Problem formulated here
class FG_eval
{
  public:
	typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
	void operator()(ADvector &fg, const ADvector &x)
	{
		assert(fg.size() == 1);
		assert(x.size() == 7);

		//data (normally would be imported, but this time hand generated)
		int data[14][4] = {
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

		//populate cost and price matrices
		//data needs to be sorted by time, then by option ID
		int prices[4][7] = {{0}}; //price matrix
		int counts[4][7] = {{0}}; //count matrix
		int t;
		int i;
		for (int s = 0; s < 14; s++)
		{
			t = data[s][0];
			i = data[s][1];
			prices[t][i] = data[s][2];
			counts[t][i] = data[s][3];
		}

		//calculating utility matrix
		AD<double> utils[4][7] = {{0}};
		for (int t = 0; t < 4; t++)
		{
			for (int i = 0; i < 7; i++)
			{
				if (counts[t][i] > 0)
				{
					utils[t][i] = x[i - 1] + x[6] * prices[t][i];
				}
			}
		}

		//calculating prob matrix
		//first need to calculate avg prob per time period
		AD<double> avg_p[4]; // array avg prob over all options in a time period
		for (int t = 0; t < 4; t++)
		{
			AD<double> sumprob = 0;
			for (int i = 0; i < 7; i++)
			{
				sumprob += exp(utils[t][i]);
			}
			avg_p[t] = sumprob;
		}
		//calculating matrix of probabilities
		AD<double> probs[4][7] = {{0}};
		for (int t = 0; t < 4; t++)
		{
			for (int i = 0; i < 7; i++)
			{
				probs[t][i] = exp(utils[t][i]) / avg_p[t];
			}
		}

		//building LL
		for (int t = 0; t < 4; t++)
		{
			for (int i = 0; i < 7; i++)
			{
				// f(x)
				fg[0] -= counts[t][i] * log(probs[t][i]);
			}
		}

		return;
	}
};
} // namespace

// m-step optimization
int optimize()
{
	bool ok;

	size_t i;
	typedef CPPAD_TESTVECTOR(double) Dvector;

	// number of independent variables (domain dimension for f and g)
	size_t nx = 7;
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
	options += "Integer print_level  0\n";
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
	// load data and preprocessing
	int **sigma_matrix = import_prefs("data/hotel_5/PrefListsBuyUpH5.csv");
	int **avail_matrix = import_availability("data/hotel_5/AvailabilityH5.csv");
	int *trans_vec = import_transactions("data/hotel_5/TransactionsH5.csv");
	int **mu_matrix = build_mu_mat(sigma_matrix, avail_matrix, trans_vec);
	// Data import debugging prints ##################################################
	{
		// printMatrix("PREFERENCE MATRIX:", sigma_matrix, 10, n_options + 1, 3);
		// printMatrix("AVAILABILITY MATRIX:", avail_matrix, 10, n_options, 3);
		// printVector("TRANSACTION VECTOR:", trans_vec, 10, 3);
		printMatrix("MU MATRIX:", mu_matrix, 20, n_types, 3);
	}

	// initialize x vector to some arbitrary values
	double x[n_types];
	std::fill_n(x, n_types, 1);

	// EM loop starts here (currently running once for testing)
	bool done = 0;
	while (!done)
	{
		// E step:
		double *m_vec = initialize_m_vec();
		double **p_sigma_matrix = build_cust_type_probs(x, mu_matrix);
		printMatrix("P_SIGMA MATRIX:", p_sigma_matrix, n_times, n_types, 5);

		done = 1;
	}

	int ok;
	ok = optimize();

	std::cout << "TEST DONE" << std::endl;
	return ok;
	//Compatability code for ipopt ##################################################
}