#include <random>
#include <chrono>
#include "common.h"
using namespace std;

// Generate an integer uniformly in range [a, b]
int Common::uniform_int_gen(int a, int b) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	uniform_int_distribution<int> distribution(a, b);
	return distribution(generator);
}

// Generate a real number uniformly in range [a, b)
double Common::uniform_real_gen(double a, double b) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	uniform_real_distribution<double> distribution(a, b);
	return distribution(generator);
}

double Common::normal_gen(double mean, double stddev) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	
	normal_distribution<double> distribution(mean, stddev);
	return distribution(generator);
}

// The DBF* function of a reservation in an interval of length t.
// See paper "The federated scheduling of constrained-deadline sporadic
// DAG task systems" by Baruah, DATE 2015.
double Common::dbf_approx(double work, double deadline, double period, double t) {
	if (t < deadline) {
		return 0;
	} else {
		double val = work + (work * (t - deadline))/period;
		return val;
	}
}
