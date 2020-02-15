#ifndef __NAIVE_GREEDY_SCHEDULER
#define __NAIVE_GREEDY_SCHEDULER

#include <queue>
#include <vector>
#include <string>
#include "dag_task.h"

struct ReadySubtask {
	uint_t id;
	uint_t ready_time;
};

// Ready subtask with earliest ready time is popped first.
struct ReadyTimeCompare {
	bool operator()(const ReadySubtask& a, const ReadySubtask& b) {
		return a.ready_time > b.ready_time;
	}
};

struct FinishedSubtask {
	uint_t id;
	uint_t start_time;
};

// Contain a list of subtasks scheduled on this processor.
struct Processor {
	Processor() {
		id = 0;
		avail_time = 0;
	}
	
	// ID of this processor.
	uint_t id;

	// Next available time of the processor.
	uint_t avail_time;
	
	// Subtasks are stored in order they are scheduled.
	list<FinishedSubtask> subtasks;
};

struct AvailTimeCompare {
	bool operator()(const Processor& a, const Processor& b) {
		return a.avail_time > b.avail_time;
	}
};

class NaiveGreedyScheduler {
public:
	NaiveGreedyScheduler(const DagTask&);
	NaiveGreedyScheduler();
	const vector<Processor>& get_schedule() const {
		return procs_;
	}
	
	void write_schedule(string) const;
	void unrestricted_schedule(const DagTask&);
	
protected:
	uint_t schedule();
	void build_trivial_schedule();
	void add_ready_subtasks(priority_queue<ReadySubtask, vector<ReadySubtask>, ReadyTimeCompare>&, vector<uint_t>&, uint_t, uint_t);
	
protected:
	DagTask task_;
	list_t depend_list_;
	uint_t num_procs_;
	vector<Processor> procs_;
};

#endif
