// Base class for the proposed federated and Baruah's federated algorithms.
#ifndef __FEDERATED_COMMON
#define __FEDERATED_COMMON

#include "taskset.h"
#include "common.h"

// Representing partition in the proposed federated scheduling.
struct ProposedPartition {
	ProposedPartition() : total_load{0} {}
	
	// Total density in case EdfSchedTest::DENSITY is used.
	// Total utilization in case EdfSchedTest::DBF_APPROX is used.
	double total_load;

	// List of light tasks assigned to this processor.
	vector<SequentialTask> tasks;
};

// For use in worst-fit heuristic.
// Partitions with smaller load are popped first.
struct WorstFitCompare {
	bool operator()(const ProposedPartition& a, const ProposedPartition& b) {
		return a.total_load > b.total_load;
	}
};

// For use in best-fit heuristic.
struct BestFitCompare {
	bool operator()(const ProposedPartition& a, const ProposedPartition& b) {
		return a.total_load < b.total_load;
	}
};

// Sort light tasks in non-decreasing order of relative deadlines.
struct NonDecreasingDeadline {
	bool operator()(const SequentialTask& a, const SequentialTask& b) {
		return a.get_deadline() < b.get_deadline();
	}
};

// Sort light tasks in non-increasing order of loads.
struct SortWorstFit {
	bool operator()(const SequentialTask& a, const SequentialTask& b) {
		return a.get_load() > b.get_load();
	}
};

// Sort light tasks in non-decreasing order of loads.
struct SortBestFit {
	bool operator()(const SequentialTask& a, const SequentialTask& b) {
		return a.get_load() < b.get_load();
	}
};

class FederatedScheduler {
public:
	FederatedScheduler(const Taskset& taskset, Heuristic h, EdfSchedTest test);
	bool is_schedulable();

protected:
	virtual bool schedule_heavy_tasks() = 0;
	bool schedule_light_tasks();
	bool worst_fit_partition();
	bool best_fit_partition();
	bool first_fit_partition();
	bool dbf_approx_partition(ProposedPartition&, const SequentialTask&);

protected:
	Taskset tset_;
	uint_t avail_procs_;
	uint_t heavy_procs_;
	Heuristic heuristic_;
	EdfSchedTest sched_test_;
};

#endif
