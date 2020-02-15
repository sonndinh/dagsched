// Implementation for the stretching algorithm in paper "A Stretching
// Algorithm for Parallel Real-Time DAG Tasks on Multiprocessor Systems"
// by Qamhieh et al., RTNS 2014.
#ifndef __STRETCHING_ALGO
#define __STRETCHING_ALGO

#include <vector>
#include "taskset.h"
#include "naive_greedy_scheduler.h"

struct Segment {
	double start; // Start time of the segment
	double len; // Length of the segment
	vector<uint_t> subtasks; // Subtasks in this segments
};

using mts_t = vector<Segment>;

// Remaining threads after stretching.
struct SequentialThread {
	double offset;
	double work;
	double deadline;
	uint_t period; // Inherit the period of DAG task.
};

// Scheduling method for sequential threads, G-EDF or P-EDF.
enum class SchedulingAlgo {GEDF, PEDF};

// Representation for partitions of P-EDF.
struct StretchingPartition {
	StretchingPartition() : total_load{0} {}

	// Total utilization of the threads assigned so far.
	double total_load;
	vector<SequentialThread> threads;
};

// The partition with smallest available space is on top.
struct StretchingBestFit {
	bool operator()(const StretchingPartition& a, const StretchingPartition& b) {
		return a.total_load < b.total_load;
	}
};

// Sort threads in non-decreasing order of relative deadlines.
struct StretchingNonDecreasingDeadline {
	bool operator()(const SequentialThread& a, const SequentialThread& b) {
		return a.deadline < b.deadline;
	}
};

// Stretch DAG tasks and apply G-EDF on the set of resulting
// constrained-deadline sequential tasks (each with intermediate
// offset and deadline).
class StretchingAlgo {
public:
	StretchingAlgo(const Taskset& taskset, SchedulingAlgo algo);
	bool is_schedulable();
	
private:
	void build_mts(const DagTask&, const NaiveGreedyScheduler&);
	void stretch(uint_t, mts_t&);
	bool schedule_threads_gedf(uint_t);
	bool schedule_threads_pedf(uint_t);
	bool dbf_approx_partition(StretchingPartition&, const SequentialThread&);
	
	// For testing. Print the MTS representation of the first heavy task.
	void print_mts() const;
	// For testing. Print the threads of the first heavy task after stretching.
	void print_stretched_threads() const;

private:
	Taskset tset_;
	SchedulingAlgo algo_;

	// Set of Multi-Threaded Segment representations, each for a heavy task.
	vector<mts_t> mtss_;
	
	// Set of sequential threads after stretching.
	// Each inner vector contains threads for a particular DAG task.
	vector<vector<SequentialThread>> threads_;
};

#endif
