#ifndef __RAND_DAG
#define __RAND_DAG

#include <vector>
#include <utility>
#include "common.h"
using namespace std;

class DagGen {	
public:
	void gen_adj_matrix(matrix_t &am, uint_t size);
	void gen_erdos_matrix(matrix_t &am, uint_t size, double p);
	uint_t gen_node_lengths(vector<uint_t> &wcets, uint_t size, uint_t min, uint_t max);
	uint_t calc_span(const matrix_t &am, const vector<uint_t> &wcets,
					 vector<pair<uint_t, int>> &dist, vector<uint_t> &crit_path);
	
	uint_t calc_subdag_work(uint_t node, const matrix_t &am, const vector<uint_t> &wcets);
	pair<vector<uint_t>, uint_t> calc_longest_path(uint_t sub_id, const matrix_t &am, const vector<uint_t> &wcets);
	void test();

private:
	bool acyclic_check(matrix_t &am, uint_t size, uint_t start, uint_t end);
	bool weak_conn_test(matrix_t &am);
	int weak_conn_add(matrix_t &am);

public:
	// Unused methods.
	void test_DAG(matrix_t &am, vector<uint_t> &wcets);
	uint_t calc_span_old(const matrix_t &am, const vector<uint_t> &wcets);
	uint_t gen_periods_NP(vector<uint_t> &wcets, uint_t size, uint_t min, uint_t max);
};

#endif
