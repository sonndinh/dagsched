#include <random>
#include <iostream>
#include <cassert>
#include "dag_task.h"

// @size: number of nodes in the DAG.
// @wcet_min: lower bound for wcet of a node.
// @wcet_max: upper bound for wcet of a node.
DagTask::DagTask(uint_t size, uint_t wcet_min, uint_t wcet_max, double p) : Task(), size_{size} {
	gen_task(size, wcet_min, wcet_max, p);
	// Generate task's deadline.
	gen_deadline();

	// For this constructor, just set period equal to the deadline.
	period_ = deadline_;
}

void DagTask::iterative_gen_task_with_util(uint_t size, uint_t wcet_min, uint_t wcet_max, double p, double util) {
	// Keep generating until we get a valid period, i.e., period > span.
	// Must clear old contents for the related data structures.
	uint_t iters = 0;
	while (true) {
		gen_task(size, wcet_min, wcet_max, p);
		period_ = (uint_t)ceil(work_/util);
		if (period_ > span_) {
			if (period_ <= work_) { // Utilization >= 1.0
				// Generate deadline uniformly in the range [span_+1, period_] for constrained task.
				deadline_ = Common::uniform_int_gen(span_+1, period_);
			} else { // Utilization < 1.0
				deadline_ = Common::uniform_int_gen(span_, work_);
			}
			break;
		}
		++iters;
		if (iters % 100 == 0) {
			cout << "Attempted generating DAG " << iters << " times ..." << endl;
		}
	}	
}

// Generate DAG task with a given utilization (and thus period).
// Note that we are generating heavy task, thus @util MUST >= 1.0.
DagTask::DagTask(uint_t size, uint_t wcet_min, uint_t wcet_max, double p, double util) : Task(), size_{size} {
	iterative_gen_task_with_util(size, wcet_min, wcet_max, p, util);
}

// Similar to the above constructor except that the number of nodes
// is uniformly chosen at random in [size_min, size_max].
DagTask::DagTask(uint_t size_min, uint_t size_max, uint_t wcet_min, uint_t wcet_max, double p, double util) : Task() {
	uint_t size = Common::uniform_int_gen(size_min, size_max);
	//	cout << "Number of generated subtasks: " << size << endl;
	size_ = size;
	iterative_gen_task_with_util(size, wcet_min, wcet_max, p, util);
}

// For generating task with a given density.
// @dummy: an additional argument to distinguish this constructor.
DagTask::DagTask(uint_t size_min, uint_t size_max, uint_t wcet_min, uint_t wcet_max, double p, double density, bool dummy) : Task() {
	// Density MUST be > 1.0.
	assert(density > 1.0);
	uint_t size = Common::uniform_int_gen(size_min, size_max);
	size_ = size;
	iterative_gen_task_with_density(size, wcet_min, wcet_max, p, density);
}

void DagTask::iterative_gen_task_with_density(uint_t size, uint_t wcet_min, uint_t wcet_max, double p, double density) {
	uint_t iters = 0;
	// Keep generating until deadline >= span.
	// Must clear old contents for the related data structures.
	while (true) {
		gen_task(size, wcet_min, wcet_max, p);
		deadline_ = (uint_t)ceil(work_/density);
		if (deadline_ >= span_) {
			// Just set period equal to deadline for now.
			// Caller MUST call change_deadline_to_period_ratio()
			// to update period to a valid value.
			period_ = deadline_;
			break;
		}
	
		++iters;
		if (iters % 100 == 0) {
			cout << "Attempted generating DAG " << iters << " times ..." << endl;
		}
	}
}

// Change the period of the task so that deadline/period falls in
// the range [min, max).
void DagTask::change_deadline_to_period_ratio(double min, double max) {
	double beta = Common::uniform_real_gen(min, max);
	period_ = (uint_t)ceil(deadline_/beta);
	//	cout << "Heavy task -- Beta " << beta << ". Deadline: " << deadline_ << ". Period: " << period_ << endl;
}

void DagTask::gen_task(uint_t size, uint_t wcet_min, uint_t wcet_max, double p) {
	// We use Erdos-Renyi method for generating DAG structure.
	DagGen gen_dag;
	gen_dag.gen_erdos_matrix(adj_matrix_, size, p);
	
	// Generate nodes's wcets.
	work_ = gen_dag.gen_node_lengths(wcets_, size_, wcet_min, wcet_max);

	// Compute span of the DAG. Note that we MUST call this
	// method after generating node lengths.
	span_ = gen_dag.calc_span(adj_matrix_, wcets_, distances_, crit_path_);

	// Compute sub-DAGs's works.
	subdag_works_.resize(size);
	for (int i = 0; i < size; ++i) {
		subdag_works_[i] = gen_dag.calc_subdag_work(i, adj_matrix_, wcets_);
	}

	// Compute longest paths.
	longest_paths_.resize(size);
	for (int i = 0; i < size; ++i) {
		longest_paths_[i] = gen_dag.calc_longest_path(i, adj_matrix_, wcets_);
	}
	
	// Construct adjacency list.
	compute_adj_list();

	// Construct dependency graph.
	compute_depend_list();
}

// Set up an object using an existing DAG. 
DagTask::DagTask(matrix_t adj_matrix, vector<uint_t> wcets, uint_t period, uint_t deadline) :
	Task(0, 0, period, deadline), size_{} {
	size_ = wcets.size();
	adj_matrix_ = adj_matrix;
	wcets_ = wcets;

	for (uint_t wcet : wcets) {
		work_ += wcet;
	}
	
	DagGen gen_dag;
	span_ = gen_dag.calc_span(adj_matrix_, wcets_, distances_, crit_path_);
	
	// Compute sub-DAGs's works.
	subdag_works_.resize(size_);
	for (int i = 0; i < size_; ++i) {
		subdag_works_[i] = gen_dag.calc_subdag_work(i, adj_matrix_, wcets_);
	}

	// Compute longest paths.
	longest_paths_.resize(size_);
	for (int i = 0; i < size_; ++i) {
		longest_paths_[i] = gen_dag.calc_longest_path(i, adj_matrix_, wcets_);
	}
	
	// Construct adjacency list.
	compute_adj_list();

	// Construct dependency graph.
	compute_depend_list();	
}

// Generate the DAG task's deadline.
void DagTask::gen_deadline() {
	// Deadline must be at least span and smaller than the task's work,
	// since we are considering heavy DAG tasks.
	// Lower-bound for the task density.
	double lowerbound = 1.0;
	double delta = 0.1; // Largest lowerbound is 1.1.
	
	// Upper-bound for the task density.
	// In normal cases, upperbound is expected to be > 1.0.
	const double upperbound = ((double)work_)/span_;

	// Special case when work is equal to span, i.e., task can only
	// execute sequentially. In this case, just set deadline equal
	// to work (and span).
	if (upperbound == 1.0) {
		deadline_ = work_;
		return;
	}
	
	// Find a valid lowerbound.
	while (lowerbound + delta >= upperbound) {
		delta = delta / 10;
	}
	lowerbound += delta;

	// Make sure the range [lowerbound, upperbound) is valid.
	// Should always satisfy.
	assert(lowerbound < upperbound);
	
	// Generate density uniformly in range (lowerbound, upperbound).
	random_device rd;
    mt19937 gen_engine(rd());
    uniform_real_distribution<double> dist(lowerbound, upperbound);
	
	double density = dist(gen_engine);
	double temp = work_/density; // Temporary value for deadline.
	deadline_ = floor(temp);
	if (deadline_ < span_) {
		deadline_ = ceil(temp); // Make sure D_i >= L_i.
	}
	//	if (upperbound < 1.1) {
	//		cout << ". Generated deadline: " << deadline_ << ". Work: " << work_ << ". Span: " << span_ << endl;
	//	}
	//	deadline_ = ceil(temp); // Should always satisfy: D_i > L_i.
	//	if (deadline_ > work_) {
	//		deadline_ = floor(temp);
	//	}
	
	// Make sure the deadline is valid.
	//	if (!(deadline_ >= span_ && deadline_ < work_)) {
	//		cout << "Work: " << work_ << ". Span: " << span_ << ". Deadline: " << deadline_ << endl;
	//	}
	assert(deadline_ >= span_ && deadline_ < work_);
}

// Compute the adjacency list from the adjacency matrix.
// Also compute the list of sinks.
void DagTask::compute_adj_list() {
	adj_list_.clear();
	adj_list_.resize(size_);

	// Go through the adjacency matrix. By the way we construct
	// the DAG, there are only edges from nodes with smaller ids to
	// nodes with larger ids.
	for (int i = 0; i < (int)size_ - 1; ++i) {
		for (int j = i + 1; j < size_; ++j) {
			if (adj_matrix_[i][j] == 1) {
				adj_list_[i].push_back(j);
			}
		}
	}

	// Compute the list of sinks.
	sinks_.clear();
	for (int i = 0; i < size_; ++i) {
		if (adj_list_[i].empty()) {
			sinks_.push_back(i);
		}
	}
}

// Compute the dependency graph. Store a list of parent subtasks
// for each subtask. This helps us know when a subtask is ready.
// Also compute the list of sources.
void DagTask::compute_depend_list() {
	depend_list_.clear();
	depend_list_.resize(size_);

	for (int i = 0; i < (int)size_ - 1; ++i) {
		for (int j = i + 1; j < size_; ++j) {
			if (adj_matrix_[i][j] == 1) {
				depend_list_[j].push_back(i);
			}
		}
	}

	// Compute the list of sources.
	sources_.clear();
	for (int i = 0; i < size_; ++i) {
		if (depend_list_[i].empty()) {
			sources_.push_back(i);
		}
	}
}
