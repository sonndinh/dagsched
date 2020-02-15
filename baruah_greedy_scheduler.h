#ifndef __BARUAH_GREEDY_SCHEDULER
#define __BARUAH_GREEDY_SCHEDULER

#include "dag_task.h"
#include "naive_greedy_scheduler.h"

class BaruahGreedyScheduler : public NaiveGreedyScheduler {
public:
	BaruahGreedyScheduler(const DagTask&);

private:
	void baruah_schedule();
	void reset();
};

#endif
