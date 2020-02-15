#ifndef __TASKSET
#define __TASKSET

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include "task.h"
#include "dag_task.h"

// Class for representing light tasks (DAG tasks with densities <= 1.0)
// and container tasks in semi-federated scheduling.
// NOTE: Do NOT call get_work(), get_span(), get_deadline(), get_period()
// for container tasks since they don't have meaningful values.
class SequentialTask : public Task {
public:
	enum class SeqTaskType {LIGHT, CONTAINER};
	
	// Constructor for creating light tasks which have densities <= 1.0.
	// The input @util must be < 1.0.
	SequentialTask(uint_t size_min, uint_t size_max, uint_t wcet_min,
				   uint_t wcet_max, double p, double util) : Task() {
		type_ = SeqTaskType::LIGHT;
		
		// Generate number of nodes.
		size_ = Common::uniform_int_gen(size_min, size_max);

		// Generate the DAG structure.
		DagGen gen_dag;
		//		matrix_t adj_matrix;
		gen_dag.gen_erdos_matrix(adj_matrix_, size_, p);
		//		vector<uint_t> wcets;
		work_ = gen_dag.gen_node_lengths(wcets_, size_, wcet_min, wcet_max);
		vector<pair<uint_t, int>> distances;
		//		vector<uint_t> crit_path;
		span_ = gen_dag.calc_span(adj_matrix_, wcets_, distances, crit_path_);

		// Compute period.
		period_ = ceil(work_/util);
		assert(period_ > work_);
		
		// Generate deadline in range [work_, period_]
		deadline_ = Common::uniform_int_gen(work_, period_);

		// For light DAG task, its load equals to its density.
		// Its load star is also equal to its load.
		load_ = (double)work_/deadline_;
		load_star_ = load_;
		splitable_ = false; // Always false for light tasks.
	}

	// Constructor for creating container task.
	// @load: The load bound of this container task.
	// @load_star: The larger portion of the load.
	SequentialTask(double load, double load_star) :	Task(), load_{load}, load_star_{load_star} {
		type_ = SeqTaskType::CONTAINER;
		splitable_ = true;
		
		// Don't care about other parameters.
		work_ = 0;
		period_ = 0;
		deadline_ = 0;
		span_ = 0;
	}

	// For generating light task with a given density.
	// Input @density MUST be <= 1.0.
	SequentialTask(uint_t size_min, uint_t size_max, uint_t wcet_min,
				   uint_t wcet_max, double p, double density, bool dummy) : Task() {
		assert(density <= 1.0);
		type_ = SeqTaskType::LIGHT;
		
		// Generate number of nodes.
		size_ = Common::uniform_int_gen(size_min, size_max);
		
		// Generate the DAG structure.
		DagGen gen_dag;
		//		matrix_t adj_matrix;
		gen_dag.gen_erdos_matrix(adj_matrix_, size_, p);
		//		vector<uint_t> wcets;
		work_ = gen_dag.gen_node_lengths(wcets_, size_, wcet_min, wcet_max);
		vector<pair<uint_t, int>> distances;
		//		vector<uint_t> crit_path;
		span_ = gen_dag.calc_span(adj_matrix_, wcets_, distances, crit_path_);
		
		// Compute deadline.
		deadline_ = (uint_t)ceil(work_/density);
		
		// Assign temporary value for period.
		// Caller MUST call to change_deadline_to_period_ratio()
		// to update period.
		period_ = deadline_;
		
		// For light DAG task, its load equals to its density.
		// Its load star is also equal to its load.
		load_ = (double)work_/deadline_;
		load_star_ = load_;
		splitable_ = false; // Always false for light tasks.
	}
	
	// Change the period of this task.
	void change_deadline_to_period_ratio(double min, double max) {
		double beta = Common::uniform_real_gen(min, max);
		period_ = (uint_t)ceil(deadline_/beta);
		//		cout << "Light task -- Beta: " << beta << ". Deadline: " << deadline_ << ". Period: " << period_ << endl;
	}
	
	enum SeqTaskType get_type() const {
		return type_;
	}

	bool is_splitable() const {
		return splitable_;
	}
	
	double get_load() const {
		return load_;
	}

	void set_load(double load) {
		load_ = load;
	}

	double get_load_star() const {
		return load_star_;
	}

	uint_t get_size() const {
		return size_;
	}

	const matrix_t& get_adj_matrix() const {
		return adj_matrix_;
	}

	const vector<uint_t>& get_node_lengths() const {
		return wcets_;
	}
	
	const vector<uint_t>& get_crit_path() const {
		return crit_path_;
	}
	
private:
	// Whether this is a light task or container task.
	SeqTaskType type_;
	// Light task cannot be split.
	// A container task in semi-federated can be split only once.
	// After that both parts cannot be split further.
	bool splitable_;
	// Update: Not really need @splitable_ attribute. Just leave it for now.
	
	// Load (also, density) of light task and load bound of
	// container task.
	double load_;
	// The larger part of the load bound of container task
	// (delta* in the paper).
	double load_star_;

	// Update (Jan 28, 2020): store the DAG structure of light tasks
	// so that we can apply Jiang et al.'s decomposition on them.
	// Number of subtasks of this DAG.
	uint_t size_;
	
	// WCETs of the subtasks.
	vector<uint_t> wcets_;

	// Adjacency matrix representation of the DAG.
	matrix_t adj_matrix_;

	// List of IDs of subtasks belonging to a critical path of the DAG.
	vector<uint_t> crit_path_;
};

class Taskset {
public:
	Taskset(uint_t m) : num_procs_{m}, total_util_{0}, heavy_util_{0} {}
	Taskset(uint_t, double, double, uint_t, string,
			uint_t, uint_t, uint_t, uint_t, double, uint_t, string);

	Taskset(uint_t, double, double, uint_t, string,
			uint_t, uint_t, uint_t, uint_t, double, uint_t, string, bool);

	void change_deadline_to_period_ratio(double min, double max);

	void add_heavy_task(const DagTask& t) {
		heavy_tasks_.push_back(t);
		total_util_ += t.get_utilization();
		heavy_util_ += t.get_utilization();
	}
	
	uint_t get_no_procs() const {
		return num_procs_;
	}

	double get_total_util() const {
		return total_util_;
	}

	double get_heavy_util() const {
		return heavy_util_;
	}

	double get_light_util() const {
		return total_util_ - heavy_util_;
	}

	const vector<DagTask>& get_heavy_tasks() const {
		return heavy_tasks_;
	}

	const vector<SequentialTask>& get_light_tasks() const {
		return light_tasks_;
	}
	
private:
	// List of heavy tasks and light tasks.
	vector<DagTask> heavy_tasks_;
	vector<SequentialTask> light_tasks_;
	
	// Number of processors in the system.
	uint_t num_procs_;
	
	// Total utilization of the task set.
	double total_util_;
	
	// Total utilization for heavy tasks.
	double heavy_util_;
	// Total utilization for light tasks equal total_util_ - heavy_util_
	// NOTE: don't need to worry about total_util_ and heavy_util_ in
	// case tasks are generated for densities. 
};

#endif
