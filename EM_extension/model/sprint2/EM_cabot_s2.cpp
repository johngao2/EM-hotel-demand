// EM implemented for Hotel 5 dataset
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm> // copy
#include <iterator>  // ostream_operator
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <boost/tokenizer.hpp>
#include <cppad/ipopt/solve.hpp>
#include "csv.h"

#define n_times 3558100 // number of time steps
#define n_options 7		// number of products
#define n_types 8		// number of customer types

// GLOBAL VARS
double m_vec[n_types];		   // m_vector, counts number of occurences of a type n arrival
double current_x_vec[n_times]; // current solution vector
double x_diff_vec[n_times];	// tracks changes in solution
double a_vec[n_times];		   // a_vector, tracks if there was an arrival in a period
double lambda;				   // arrival parameter
double alpha = 0.1;			   // regularization hyperparameter
int n_purch;				   // tracks number total number of purchases

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
void printVector(const char *text, T mat, std::size_t N, int width, int precision = 2)
{
	std::cout << text << std::endl;

	for (int i = 0; i < N; i++)
	{
		std::cout << std::setw(width) << std::fixed << std::setprecision(precision)
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
	// importing csv
	using namespace boost;
	std::string data(pref_filename);
	std::ifstream in(data.c_str());

	// error check
	if (!in.is_open())
	{
		std::cout << "ERROR: CSV NOT FOUND" << std::endl;
		std::exit(1);
	}

	// declare parsing vars
	typedef tokenizer<escaped_list_separator<char>> Tokenizer;
	std::vector<std::string> vec;
	std::string line;
	int header_tf = 1;	// helper vars to ignore header and index values
	int index_tf = 1;
	int col_counter = 0; // col counter for to select which col to insert value
	int row_counter = 0;

	// vars to store data
	int **pref_matrix = 0; // initialize
	pref_matrix = new int *[n_types];

	std::cout << "Loading preferences" << std::endl;
	// read line by line
	while (getline(in, line))
	{
		// declare stuff
		index_tf = 1;
		col_counter = 0;
		Tokenizer tok(line);
		pref_matrix[row_counter] = new int[n_options + 1]; // new row

		// ignore first line
		if (header_tf == 1)
		{
			header_tf = 0;
			continue;
		}

		// iterate over row tokens
		for (Tokenizer::iterator it(tok.begin()),
								end(tok.end());
			it != end; ++it)
		{
			// skip first value since it's index
			if(index_tf == 1){
				index_tf = 0;
				continue;
			}
			pref_matrix[row_counter][col_counter] = std::stoi(*it);
			col_counter ++;
		}
		row_counter ++;

		// printing progress
		if(row_counter % 100000 == 0){
			std::cout << row_counter << std::endl;
		}
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
	std::cout << "Done preferences" << std::endl;
	return sigma_matrix;
}

// import availability matrix
// returns avail_matrix: a n_times x n_options matrix
int **import_availability(const char *avail_filename)
{
	std::cout << "Loading avail" << std::endl;
	// importing csv
	using namespace boost;
	std::string data(avail_filename);
	std::ifstream in(data.c_str());

	// error check
	if (!in.is_open())
	{
		std::cout << "ERROR: CSV NOT FOUND" << std::endl;
		std::exit(1);
	}

	// declare parsing vars
	typedef tokenizer<escaped_list_separator<char>> Tokenizer;
	std::vector<std::string> vec;
	std::string line;
	int header_tf = 1;	// helper vars to ignore header and index values
	int index_tf = 1;
	int col_counter = 0; // col counter for to select which col to insert value
	int row_counter = 0;

	// vars to store data
	int **avail_matrix = 0; // initialize
	avail_matrix = new int *[n_times];

	// read line by line
	while (getline(in, line))
	{
		// declare stuff
		index_tf = 1;
		col_counter = 0;
		Tokenizer tok(line);
		avail_matrix[row_counter] = new int[n_options]; // new row

		// ignore first line
		if (header_tf == 1)
		{
			header_tf = 0;
			continue;
		}

		// iterate over row tokens
		for (Tokenizer::iterator it(tok.begin()),
								end(tok.end());
			it != end; ++it)
		{
			// skip first value since it's index
			if(index_tf == 1){
				index_tf = 0;
				continue;
			}
			avail_matrix[row_counter][col_counter] = std::stoi(*it);
			col_counter ++;
		}
		row_counter ++;

		// printing progress
		if(row_counter % 100000 == 0){
			std::cout << row_counter << std::endl;
		}
	}
	std::cout << "Done avail" << std::endl;
	return avail_matrix;
}

// import transaction vector
// returns trans_vector: a length T vector
int *import_transactions(const char *trans_filename)
{
	std::cout << "Loading trans" << std::endl;
	// importing csv
	using namespace boost;
	std::string data(trans_filename);
	std::ifstream in(data.c_str());

	// error check
	if (!in.is_open())
	{
		std::cout << "ERROR: CSV NOT FOUND" << std::endl;
		std::exit(1);
	}

	// declare parsing vars
	typedef tokenizer<escaped_list_separator<char>> Tokenizer;
	std::vector<std::string> vec;
	std::string line;
	int header_tf = 1;	// helper vars to ignore header and index values
	int index_tf = 1;
	int row_counter = 0;

	// vars to store data
	int *trans_vec = new int[n_times];

	// read line by line
	while (getline(in, line))
	{
		// declare stuff
		index_tf = 1;
		Tokenizer tok(line);

		// ignore first line
		if (header_tf == 1)
		{
			header_tf = 0;
			continue;
		}

		// iterate over row tokens
		for (Tokenizer::iterator it(tok.begin()),
								end(tok.end());
			it != end; ++it)
		{
			// skip first value since it's index
			if(index_tf == 1){
				index_tf = 0;
				continue;
			}
			trans_vec[row_counter] = std::stoi(*it);
		}
		row_counter ++;

		// printing progress
		if(row_counter % 100000 == 0){
			std::cout << row_counter << std::endl;
		}
	}
	std::cout << "Done trans" << std::endl;
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
				bool available = avail_matrix[t][k - 1]; // need to subtract 1 from k due to 1-indexed sigma matrix
				if (available == 1 && k != j_t)			 // if item is available and not selected
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

// literally just counts number of purchases
void count_purchases(int *trans_vec)
{
	for (int t = 0; t < n_times; t++)
	{
		if (trans_vec[t] != 0)
		{
			n_purch++;
		}
	}
}

// build p_sigma matrix (calculates customer type probabilities
// based on the data and previous guess for xi)
// p_sigma is confusingly denoted as x_it in the pseudocode
// resurns p_sigma matrix: n_times x n_types matrix
double **build_cust_type_probs(int **mu_matrix)
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
			sumprob += mu_matrix[t][i] * current_x_vec[i];
		}

		// update type probability if that type is compatible
		for (int i = 0; i < n_types; i++)
		{
			if (mu_matrix[t][i] == 1) // set was compatible, so update this value
			{
				p_sigma_matrix[t][i] = current_x_vec[i] / sumprob;
			}
			else
			{
				p_sigma_matrix[t][i] = 0;
			}
		}
	}
	return p_sigma_matrix;
}

// updates a_t estimates based on compatible types and transaction data
void update_arrival_estimates(int *trans_vec, int **mu_matrix)
{
	for (int t = 0; t < n_times; t++)
	{
		if (trans_vec[t] != 0) // case 1: definitely arrival, item was purchased
		{
			a_vec[t] = 1;
		}
		else
		{
			int num_compat_types = 0;
			for (int i = 0; i < n_types; i++)
			{
				num_compat_types += mu_matrix[t][i];
			}
			if (num_compat_types == 0) // case 2: definitely non-arrival, no compatible types
			{
				a_vec[t] = 0;
			}
			else // case 3: might have been an arrival, update a_t with probability
			{
				double sum_compat_probs = 0; // sum of probabilities of compatible types
				for (int i = 0; i < n_types; i++)
				{
					sum_compat_probs += current_x_vec[i] * mu_matrix[t][i];
				}
				a_vec[t] = (lambda * sum_compat_probs) / (lambda * sum_compat_probs + (1 - lambda));
			}
		}
	}
}

// estimate m vector, last part of e-step
// returns m_vec: length N
void estimate_m_vec(double **p_sigma_matrix)
{
	for (int i = 0; i < n_types; i++)
	{
		m_vec[i] = 0;
		for (int t = 0; t < n_times; t++)
		{
			m_vec[i] += a_vec[t] * p_sigma_matrix[t][i];
		}
	}
}

// Problem formulated here
class FG_eval
{
public:
	typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
	void operator()(ADvector &fg, const ADvector &x)
	{
		assert(fg.size() == 2); // 1 LL formula, 1 sum of x constraint, 2*n_types regularization constraints
		assert(x.size() == n_types + 1);

		// printVector("m_vec:", m_vec, n_types, 4, 5);

		//x[0:n_types) represents type probs, x[n_types] is lambda

		// add type probs to objective
		for (int i = 0; i < n_types; i++)
		{
			fg[0] += m_vec[i] * log(x[i]);
		}
		// add lambda terms to objective
		double sum_ambiguous_a = 0;
		for (int t = 0; t < n_times; t++)
		{
			if (a_vec[t] != 1)
			{
				sum_ambiguous_a += a_vec[t];
			}
		}

		AD<double> reg_sum = 0;
		// sum regularization terms
		for (int i = 0; i < n_types; i++)
		{
			reg_sum += pow(x[i], 2);
		}
		reg_sum = reg_sum * alpha;

		// main function to optimize
		fg[0] += (n_purch + sum_ambiguous_a) * log(x[n_types]) + ((n_times - n_purch) - sum_ambiguous_a) * log(1 - x[n_types]) - reg_sum;

		// constraint that sum of x's = 1
		for (int i = 0; i < n_types; i++)
		{
			fg[1] += x[i];
		}

		return;
	}
};
} // namespace

// m-step optimization
void m_step()
{
	bool ok;

	size_t i;
	typedef CPPAD_TESTVECTOR(double) Dvector;

	// number of independent variables n_types + 1 extra for lambda + regularizers vars
	size_t nx = n_types + 1;
	// number of constraints (range dimension for g)
	size_t ng = 1;
	// initial value of the independent variables

	Dvector xi(nx);
	for (i = 0; i < n_types + 1; i++)
	{
		xi[i] = 0.1;
	}

	// lower and upper limits for x
	Dvector xl(nx), xu(nx);
	for (i = 0; i < nx; i++)
	{
		xl[i] = 0;
		xu[i] = 1;
	}

	// lower and upper limits for g
	Dvector gl(ng), gu(ng);
	for (i = 0; i < ng; i++)
	{
		gl[i] = 0;
		gu[i] = 1;
	}

	// object that computes objective and constraints
	FG_eval fg_eval;

	// options
	std::string options;
	// turn off any printing
	options += "Integer print_level  0\n";
	options += "String  sb           yes\n";
	// scaling to maximize instead of minimize
	options += "Numeric obj_scaling_factor   -1\n";
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

	// printVector("PAST SOLUTION", current_x_vec, n_types, 3, 5);
	// printVector("NEW SOLUTION", solution.x, n_types, 3, 5);

	// update current optimal x vector and difference vector, as well as lambda
	for (int i = 0; i < n_types; i++)
	{
		x_diff_vec[i] = std::abs(current_x_vec[i] - solution.x[i]);
		current_x_vec[i] = solution.x[i];
	}
	lambda = solution.x[n_types];

	std::cout << "CURRENT OBJECTIVE VALUE: " << solution.obj_value << std::endl;
	// std::cout << "CURR LAMBDA: " << solution.x[n_types] << std::endl;

	// printVector("DIFFERENCE", x_diff_vec, n_types, 3, 5);
}

// alternate m_step function for testing
void closed_form_m_step()
{
	// calculate sum of m vec
	double m_sum;
	for (int i = 0; i < n_types; i++)
	{
		m_sum += m_vec[i];
	}

	// update x vec using closed form solution
	double new_x;
	for (int i = 0; i < n_types; i++)
	{
		new_x = m_vec[i] / m_sum;
		x_diff_vec[i] = std::abs(new_x - current_x_vec[i]);
		current_x_vec[i] = new_x;
	}

	// calculate objective value
	double LL;
	for (int i = 0; i < n_types; i++)
	{
		LL += m_vec[i] * log(current_x_vec[i]);
	}
}

// calculates different LL function using equation (2)
double real_LL(int **mu_matrix)
{
	double LL = 0;
	for (int t = 0; t < n_times; t++)
	{
		double temp = 0; // temp variable to store xi's for one time period
		for (int i = 0; i < n_types; i++)
		{
			temp += current_x_vec[i] * mu_matrix[t][i];
		}
		temp = log(temp);
		LL += temp;
	}

	return LL;
}

int main()
{
	// load data and preprocessing
	int **sigma_matrix = import_prefs("../../../data/cabot_data/sprint_1/types_sprint1.csv");
	int **avail_matrix = import_availability("../../../data/cabot_data/sprint_1/availability_sprint1.csv");
	int *trans_vec = import_transactions("../../../data/cabot_data/sprint_1/transactions_pre_sprint1.csv");
	int **mu_matrix = build_mu_mat(sigma_matrix, avail_matrix, trans_vec);

	// Data import debugging prints ##################################################
	{
		printMatrix("PREFERENCE MATRIX:", sigma_matrix, 8, n_options + 1, 3);
		printMatrix("AVAILABILITY MATRIX:", avail_matrix, 10, n_options, 3);
		printVector("TRANSACTION VECTOR:", trans_vec, 10, 3);
		printMatrix("MU MATRIX:", mu_matrix, 10, n_types, 3);
	}

	// init: set a_vec to 0 x_vec to 1/N, lambda to 0.5, count purchases
	std::fill_n(a_vec, n_times, 0);
	std::fill_n(current_x_vec, n_types, 1.0 / n_types);
	lambda = 0.3;
	count_purchases(trans_vec);
	std::cout << "NUM_PURCHASES: " << n_purch << std::endl;
	// initialization debugging prints:
	{
		// printVector("A_VEC", a_vec, 10, 5);
		// printVector("current_x_vec", current_x_vec, 10, 5);
	}
	double maxdiff;

	// EM loop starts here
	bool done = 0;
	while (!done)
	{
		// E step:
		// update cust type probs
		double **p_sigma_matrix = build_cust_type_probs(mu_matrix);
		// update a_t predictions
		update_arrival_estimates(trans_vec, mu_matrix);
		//printVector("A_VEC", a_vec, n_times, 5);
		estimate_m_vec(p_sigma_matrix);

		// M step:
		m_step();

		// find max difference of solution, exit loop if small enough
		maxdiff = *std::max_element(x_diff_vec, x_diff_vec + n_types);
		if (1)
		{
			done = 1;
		}
		// printVector("CURRENT SOLUTION", current_x_vec, n_types, 3, 5);
	}
	printVector("FINAL X_VEC", current_x_vec, n_types, 5, 5);
	std::cout << "FINAL LAMBDA: " << lambda << std::endl;
	std::cout << "TEST DONE" << std::endl;
	return 1;
}