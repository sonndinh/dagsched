#include <cmath>
#include <iostream>
#include "baruah_greedy_scheduler.h"

BaruahGreedyScheduler::BaruahGreedyScheduler(const DagTask& t) : NaiveGreedyScheduler(t) {
	uint_t work = task_.get_work();
	uint_t deadline = task_.get_deadline();

	// Start with a minimum number of processors.
	num_procs_ = ceil((double)work/deadline);
	reset();

	// Schedule the task.
	baruah_schedule();
}

// Find a smallest number of processors that can schedule the task.
void BaruahGreedyScheduler::baruah_schedule() {
	uint_t deadline = task_.get_deadline();

	// Keep increasing the number of processors until schedulable.
	while (true) {
		uint_t finish_time = schedule();
		// Successfully schedule the task.
		if (finish_time <= deadline)
			return;
		
		// Try adding 1 more processor.
		++num_procs_;
		reset();
	}
}

// Reset the data for the processors assigned to this task.
void BaruahGreedyScheduler::reset() {
	depend_list_ = task_.get_depend_list();
	procs_.resize(num_procs_);
	
	for (int i = 0; i < num_procs_; ++i) {
		procs_[i].id = i;
		procs_[i].avail_time = 0;
		procs_[i].subtasks.clear();
	}
}
