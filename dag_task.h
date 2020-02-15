#ifndef __DAG_TASK
#define __DAG_TASK

#include "task.h"
#include "rand_dag.h"
using namespace std;

// Implementation for heavy DAG tasks, which have densities > 1.0.
// We implement light DAG tasks (densities <= 1.0) in class SequentialTask.
class DagTask : public Task {
public:
	DagTask() {}
	DagTask(uint_t size, uint_t wcet_min, uint_t wcet_max, double p);
	DagTask(matrix_t adj_matrix, vector<uint_t> wcets, uint_t period, uint_t deadline);
	DagTask(uint_t size, uint_t wcet_min, uint_t wcet_max, double p, double util);
	DagTask(uint_t size_min, uint_t size_max, uint_t wcet_min, uint_t wcet_max, double p, double util);
	DagTask(uint_t size_min, uint_t size_max, uint_t wcet_min, uint_t wcet_max, double p, double density, bool dummy);

	void change_deadline_to_period_ratio(double min, double max);
	
	uint_t get_size() const {
		return size_;
	}
	
	const matrix_t& get_adj_matrix() const {
		return adj_matrix_;
	}

	const list_t& get_adj_list() const {
		return adj_list_;
	}

	const list_t& get_depend_list() const {
		return depend_list_;
	}

	const vector<uint_t>& get_node_lengths() const {
		return wcets_;
	}
	
	const vector<uint_t>& get_crit_path() const {
		return crit_path_;
	}

	const vector<uint_t>& get_subdag_works() const {
		return subdag_works_;
	}

	const vector<uint_t>& get_longest_paths() {
		if (path_lengths_.empty()) {
			// Only compute the path lengths array on the first call.
			path_lengths_.resize(size_);
			for (uint_t i = 0; i < size_; ++i) {
				uint_t length = longest_paths_[i].second;
				path_lengths_[i] = length;
			}
		}
		
		return path_lengths_;
	}
	
	const vector<uint_t>& get_sources() const {
		return sources_;
	}

	const vector<uint_t>& get_sinks() const {
		return sinks_;
	}
	
private:
	void gen_task(uint_t, uint_t, uint_t, double);
	void iterative_gen_task_with_util(uint_t, uint_t, uint_t, double, double);
	void iterative_gen_task_with_density(uint_t, uint_t, uint_t, double, double);
	void compute_adj_list();
	void compute_depend_list();
	void gen_deadline();
	
private:
	// DAG structure parameters.
	// Subtasks are numbered from 0 to size - 1.
	uint_t size_;

	// WCETs of the nodes.
	vector<uint_t> wcets_;

	// Adjacency matrix representation of the DAG.
	matrix_t adj_matrix_;
	
	// Adjacency list representation of the DAG.
	list_t adj_list_;

	// Dependency list for the DAG.
	list_t depend_list_;

	// List of IDs of sources.
	vector<uint_t> sources_;

	// List of IDs of sinks.
	vector<uint_t> sinks_;
	
	// Longest distances to individual nodes from a source.
	// The distance to a node includes the wcet of the node itself.
	// Node IDs (starting from 0) are used as indexes for their data.
	// Each pair contains the longest distance to the node as the 1st value,
	// and ID of its parent in the longest path as the 2nd value.
	// In case the node has no parent (i.e., a source node), the 2nd value is -1.
	vector<pair<uint_t, int>> distances_;

	// A critical-path of the DAG.
	// The vector contains IDs of nodes belonging to the critical-path in order. 
	// In case there are more than one critcal path exists, we pick an arbitrary critical path. 
	vector<uint_t> crit_path_;

	// Each value is the total work of a sub-DAG rooted at a particular node.
	// The IDs of the nodes are used as the indexes on this vector, e.g.,
	// the value for sub-DAG rooted at node 0 is the 1st entry, and so on.
	vector<uint_t> subdag_works_;

	// Store a longest path from each subtask to a sink subtask.
	// Use the subtask IDs to access their longest paths.
	// Each longest path contains the IDs of subtasks in that path
	// in its order and the length of the path.
	vector<pair<vector<uint_t>, uint_t>> longest_paths_;
	vector<uint_t> path_lengths_; // Only contain lengths of those paths.
};

#endif
