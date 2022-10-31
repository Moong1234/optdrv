#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

// Returns European payoff at the maturity date
double payoff(double price, double strike, bool put = true)
{
    return std::max(put * (strike - price), abs(put - 1) * (price - strike));
}

// Pricing and hedging an American option.
//
// INPUT: input.txt
// Arguments should be given as a single line in input.txt
// "S_0,T,K,r,sigma,n_steps"
// double S_0: initial stock price
// double T: maturity date
// double K: strike price
// double r: risk free rate
// double sigma: volatility
// int N : number of steps in your binomial tree
//
// RETURN X_0, Delta_0;
// X_0: value of the option at time 0
// Delta_0: number of shares of S for replicating portfolio at time 0
int main(int argc, char *argv[])
{
    // Read input data
    std::fstream myfile("input.txt", std::ios_base::in);
    std::string tmp;
    std::vector<std::string> v;
    while (myfile.good())
    {
        std::getline(myfile, tmp, ',');
        v.push_back(tmp);
    }

    const double S = std::stod(v[0]);
    const double T = std::stod(v[1]);
    const double K = std::stod(v[2]);
    const double R = std::stod(v[3]);
    const double SD = std::stod(v[4]);
    const int N = std::stoi(v[5]);

    // const for calculation
    const double delta_t = T / N;
    const double u = exp(SD * sqrt(delta_t));
    const double d = exp(-SD * sqrt(delta_t));
    const double p = (exp(R * delta_t) - d) / (u - d);
    const double q = 1 - p;

    std::vector<double> curr_value(N, 0);
    std::vector<double> next_value(N + 1, 0);

    double discount, hold, exec, curr_price;

    for (int i = 0; i < N + 1; i++)
        next_value[i] = payoff(pow(u, i) * pow(d, (N - i)) * S, K);

    for (int n = N; n > 0; n--)
    {
        for (int i = 0; i < n; i++)
        {
            discount = exp(-R * delta_t);
            curr_price = pow(u, i) * pow(d, ((n - 1) - i)) * S;
            exec = payoff(curr_price, K);
            hold = discount * (p * next_value[i + 1] + q * next_value[i]);
            curr_value[i] = std::max(exec, hold);
        }

        if (n > 1)
        {
            next_value = curr_value;
            curr_value.pop_back();
        }
    }

    double Delta_zero = (next_value[1] - next_value[0]) / (u * S - d * S);

    std::cout << curr_value[0] << ',' << Delta_zero << std::endl;

    return 0;
}
