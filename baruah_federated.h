#ifndef __BARUAH_FEDERATED
#define __BARUAH_FEDERATED

#include "federated_common.h"

class BaruahFederated : public FederatedScheduler {
public:
	BaruahFederated(const Taskset& taskset, Heuristic h, EdfSchedTest test);
	
private:
	bool schedule_heavy_tasks();
};

#endif
