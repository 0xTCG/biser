
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sstream>
#include <unordered_map>
#include <iostream>

using namespace std;
extern "C" double tau_c(double edit_error, int kmer_size, double MAX_ERROR,double MAX_EDIT_ERROR)
{
	// cout << MAX_EDIT_ERROR << "\t" << (MAX_ERROR - MAX_EDIT_ERROR) << endl;
	const double ERROR_RATIO = (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR;
	
	double gap_error = std::min(1.0, ERROR_RATIO * edit_error);
	double a = (1 - gap_error) / (1 + gap_error);
	double b = 1 / (2 * std::exp(kmer_size * edit_error) - 1);
	return a * b;
}

extern "C" double solve_inverse_jaccard_c(int j, int kmer_size, double MAX_ERROR,double MAX_EDIT_ERROR)
{
	if (j == 0)
	{
		return 1;
	}
	if (j == 1)
	{
		return 0;
	}
	return boost::math::tools::newton_raphson_iterate([j, kmer_size,MAX_ERROR,MAX_EDIT_ERROR](double d){
		const double ERROR_RATIO = (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR;
		double E = exp(d * kmer_size);
		return make_tuple(
			((1 - d * ERROR_RATIO) / (1 + d * ERROR_RATIO)) * (1.0 / (2 * E - 1)) - j,
			2 * (- kmer_size * E + ERROR_RATIO - 2 * ERROR_RATIO * E + E * kmer_size * pow(d * ERROR_RATIO, 2)) /
				pow((2 * E - 1) * (1 + d * ERROR_RATIO), 2)
		);
	}, 0.10, 0.0, 1.0, numeric_limits<double>::digits);
}


extern "C" double relaxed_jaccard_estimate_c(int s, int kmer_size, double  MAX_EDIT_ERROR, double MAX_ERROR, double result)
{
	
	using namespace boost::math;
	const double CI = 0.75;
	const double Q2 = (1.0 - CI) / 2; // one side interval probability
	// cout << "tau" << endl;
	result = ceil(s * tau_c(MAX_EDIT_ERROR, kmer_size,MAX_ERROR,MAX_EDIT_ERROR));
	// cout << result << endl;
	for (; result >= 0; result--) {
		
		double d = solve_inverse_jaccard_c(result / s, kmer_size,MAX_ERROR,MAX_EDIT_ERROR); // returns edit error
		
		double x = quantile(complement(binomial(s, tau_c(d, kmer_size,MAX_ERROR,MAX_EDIT_ERROR)), Q2)); // inverse binomial 
		double low_d = solve_inverse_jaccard_c(x / s, kmer_size,MAX_ERROR,MAX_EDIT_ERROR);
		if (100 * (1 - low_d) < MAX_EDIT_ERROR) {
			result++; 
			break;
		}
	}
	
	result = max(result, 0.0);
	return result;
}


extern "C" int fun1(int a){
	return a*4;
}


