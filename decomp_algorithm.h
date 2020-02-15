// Implementation of Jiang et al.'s paper "On the Decomposition-based Global
// EDF Scheduling of Parallel Real-Time Tasks", RTSS 2016.
#ifndef __DECOMP_ALGO
#define __DECOMP_ALGO

#include <list>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "taskset.h"

// Lifetime window for each node (or part of a node).
struct NodeWindow {
	uint_t id; // ID of this node
	uint_t len; // Execution time of this node.
	uint_t ready; // Earliest ready time
	uint_t finish; // Latest finish time
};

struct IncreasingFinish {
	bool operator()(const NodeWindow& a, const NodeWindow& b) {
		return a.finish > b.finish;
	}
};

struct DcomSegment {
	DcomSegment() : start{0}, end{0}, work{0} {}
	
	uint_t start; // Start time
	uint_t end; // End time
	// List of nodes (or parts of nodes) assigned to this segment.
	unordered_map<uint_t, uint_t> id_to_len;
	uint_t work; // Total work of this segment.
};

// Timing diagram for each DAG task.
struct TimingDiagram {
	// Store window of nodes in their index order.
	vector<NodeWindow> subtasks;
	double omega; // Structure characteristic value of this task.
	
	// Segments are separated by ready times of nodes.
	// The last segment ends at the end of span.
	vector<DcomSegment> segments;
};

enum class TaskType {HEAVY, LIGHT};

class DecompAlgo {
public:
	DecompAlgo(const Taskset& taskset);
	bool is_schedulable();

private:
	void decompose(const Task *);
	uint_t compute_timing_diagram(const Task*);
	list<uint_t> topological_sort(const matrix_t&);
	void dfs_visit(const matrix_t&, uint_t, list<uint_t>&, unordered_set<uint_t>&);

	void print_timing_diagram(uint_t);
	void print_segments(uint_t);

private:
	Taskset tset_;

	// Set of timing diagrams, each for a DAG task.
	vector<TimingDiagram> diagrams_;
};

#endif
