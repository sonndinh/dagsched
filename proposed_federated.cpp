#include <iostream>
#include "proposed_federated.h"

ProposedFederated::ProposedFederated(const Taskset& taskset, Heuristic h, EdfSchedTest test) : FederatedScheduler(taskset, h, test) {
}

// Call after is_schedulable() is called.
uint_t ProposedFederated::procs_to_heavy_tasks() const {
	return heavy_procs_;
}

// Schedule each heavy task using the proposed algorithm.
bool ProposedFederated::schedule_heavy_tasks() {
	const vector<DagTask>& heavy_tasks = tset_.get_heavy_tasks();
	for (int i = 0; i < heavy_tasks.size(); ++i) {
		ProposedScheduler sched(heavy_tasks[i]);
		uint_t required_cores = sched.get_num_procs();
		if (required_cores > avail_procs_)
			return false;

		avail_procs_ -= required_cores;
		heavy_procs_ += required_cores;
	}
	
	return true;
}
