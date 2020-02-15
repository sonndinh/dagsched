// Implementation for the algorithm in paper "Semi-Federated Scheduling of
// Parallel Real-Time Tasks on Multiprocessors" by Jiang et al., published at RTSS 2017.
#ifndef __SEMI_FED
#define __SEMI_FED

#include "taskset.h"

struct NonIncreasingLoadStar {
	// Used to sort delta* of sequential tasks in non-increasing order.
	bool operator()(const SequentialTask& a, const SequentialTask& b) {
		return a.get_load_star() > b.get_load_star();
	}
};

struct NonIncreasingLoad {
	// Used to sort the remaining load of container tasks after scraping.
	bool operator()(const SequentialTask &a, const SequentialTask &b) {
		return a.get_load() > b.get_load();
	}
};

// Representing a partition, i.e., processor, when doing partitioning
// in semi-federated scheduling.
struct Partition {
	Partition() {
		total_load = 0;
		total_load_star = 0;
	}

	// Sum of delta.
	double total_load;
	
	// Sum of delta*.
	double total_load_star;
	vector<SequentialTask> tasks;
};

struct SmallerLoadStar {
	bool operator()(const Partition& a, const Partition& b) {
		return a.total_load_star > b.total_load_star;
	}
};

struct SmallerLoad {
	bool operator()(const Partition& a, const Partition& b) {
		return a.total_load > b.total_load;
	}
};

class SemiFederated {
public:
	SemiFederated(const Taskset& taskset);
	bool is_schedulable();
	uint_t procs_to_heavy_tasks() const;

private:
	bool schedule_heavy_tasks();
	bool schedule_seq_tasks();
	vector<SequentialTask> scrape(Partition& par);
	
private:
	Taskset tset_;
	vector<SequentialTask> container_tasks_;
	// Total available processors. It is updated while the algorithm runs.
	uint_t avail_procs_;

	// Number of processors dedicatedly allocated to heavy tasks.
	uint_t heavy_procs_;
};

#endif
