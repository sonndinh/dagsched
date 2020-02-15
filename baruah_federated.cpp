#include "baruah_federated.h"
#include "baruah_greedy_scheduler.h"

BaruahFederated::BaruahFederated(const Taskset& taskset, Heuristic h, EdfSchedTest test) : FederatedScheduler(taskset, h, test) {
}

bool BaruahFederated::schedule_heavy_tasks() {
	const vector<DagTask>& heavy_tasks = tset_.get_heavy_tasks();
	for (int i = 0; i < heavy_tasks.size(); ++i) {
		BaruahGreedyScheduler sched(heavy_tasks[i]);
		uint_t required_cores = sched.get_schedule().size();
		if (required_cores > avail_procs_)
			return false;

		avail_procs_ -= required_cores;
		heavy_procs_ += required_cores;
	}
	
	return true;
}
