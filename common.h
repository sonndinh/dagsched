#ifndef __COMMON
#define __COMMON

#include <vector>
#include <list>

using uint_t = unsigned int;
using matrix_t = std::vector<std::vector<uint_t>>;
using list_t = std::vector<std::list<uint_t>>;

class Common {
public:
	static int uniform_int_gen(int a, int b);
	static double uniform_real_gen(double a, double b);
	static double normal_gen(double mean, double stddev);
	static double dbf_approx(double work, double deadline, double period, double t);
};

// Uniprocessor schedulability test for EDF.
// Currently support 2 methods: DBF* method and total density sufficient test.
enum class EdfSchedTest {DBF_APPROX, DENSITY};

// Heuristics for partitioning.
enum class Heuristic {WORST_FIT, BEST_FIT, FIRST_FIT};

#endif
