// Actually a logit exercise
// Works with example solver
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <math.h>

#define n_times 245
#define n_options 8
#define n_types 15

template<typename T>
void printMatrix(T mat, std::size_t N, std::size_t M) {
    cout<<"\n Printing Matrix : \n";
    for(int i = 0 ; i < N ; ++i) {
        for(int j = 0 ; j < N; ++j)
            cout<< *(*(mat+i)+j)<<" ";
        cout<<endl;
    }
    cout<<endl;
}

int main(){
	//data (normally would be imported, but this time hand generated)
		int data[n_samples][4] = {
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
		AD<double> prices[n_times][n_options] = {{0}}; //price matrix
		AD<double> counts[n_times][n_options] = {{0}}; //count matrix
		int t;
		int i;
		for (int s = 0; s < n_samples; s++)
		{
			t = data[s][0];
			i = data[s][1];
			prices[t][i] = data[s][2];
			counts[t][i] = data[s][3];
		}
}
