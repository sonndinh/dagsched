#ifndef __PROPOSED_TASKSET_SCHEDULER
#define __PROPOSED_TASKSET_SCHEDULER

#include "proposed_scheduler.h"
#include "federated_common.h"
#include "taskset.h"
#include "common.h"


// Schedule task set using the proposed algorithm for heavy tasks.
class ProposedFederated : public FederatedScheduler {
public:
	ProposedFederated(const Taskset& taskset, Heuristic h, EdfSchedTest test);
	uint_t procs_to_heavy_tasks() const;

private:
	bool schedule_heavy_tasks();
};

#endif
