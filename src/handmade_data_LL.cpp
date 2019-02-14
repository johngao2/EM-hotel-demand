#include <iostream>
#include <math.h>

#define n_times 4
#define n_options 7
#define n_samples 14

void print_data(int arr[n_samples][4])
{
    printf("ARRAY:\n");
    for (int i = 0; i < n_samples; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
}

void print_double_matrix(double arr[n_times][n_options])
{
    printf("ARRAY:\n");
    for (int i = 0; i < n_times; ++i)
    {
        for (int j = 0; j < n_options; ++j)
        {
            printf("%f ", arr[i][j]);
        }
        printf("\n");
    }
}

void print_int_matrix(int arr[n_times][n_options])
{
    printf("ARRAY:\n");
    for (int i = 0; i < n_times; ++i)
    {
        for (int j = 0; j < n_options; ++j)
        {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
}

int main(void)
{
    //fg.size() == 1
    //x.size()  == 7
    //temp decision variables for testing
    double x1 = 4; //a1
    double x2 = 1; //a2
    double x3 = 2; //a3
    double x4 = 4; //a4
    double x5 = 3; //a5
    double x6 = 7; //a6
    double x7 = 2; //b
    double theta[n_options + 1] = {x1, x2, x3, x4, x5, x6, x7};

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
    int prices[n_times][n_options] = {{0}}; //price matrix
    int counts[n_times][n_options] = {{0}}; //count matrix
    int t;
    int i;
    for (int s = 0; s < n_samples; s++)
    {
        t = data[s][0];
        i = data[s][1];
        prices[t][i] = data[s][2];
        counts[t][i] = data[s][3];
    }

    //calculating utility matrix
    double utils[n_times][n_options] = {{0}};
    int a;
    int b;
    for (int t = 0; t < n_times; t++)
    {
        for (int i = 0; i < n_options; i++)
        {
            if (counts[t][i] > 0)
            {
                a = theta[i - 1];
                b = theta[6];
                utils[t][i] = a + b * prices[t][i];
            }
        }
    }

    //calculating prob matrix
    //first need to calculate avg prob per time period
    double avg_p[n_times]; // array avg prob over all options in a time period
    for (int t=0; t<n_times; t++){
        double sumprob = 0;
        for(int i=0; i<n_options; i++){
            sumprob += exp(utils[t][i]);
        }
        avg_p[t] = sumprob;
    }
    //calculating matrix of probabilities
    double probs[n_times][n_options] = {{0}};
    for (int t=0; t<n_times; t++){
        for(int i=0; i<n_options; i++){
            probs[t][i] = exp(utils[t][i]) / avg_p[t];
        }
    }

    //building LL
    double LL;
    for (int t = 0; t < n_times; t++)
    {
        for (int i = 0; i < n_options; i++)
        {
            LL += counts[t][i] * log(probs[t][i]);
        }
    }
    printf("Current LL: %f", LL);
}
