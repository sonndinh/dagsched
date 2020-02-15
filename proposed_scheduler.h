#ifndef __PROPOSED_SCHEDULER
#define __PROPOSED_SCHEDULER

#include <vector>
#include <list>
#include <queue>
#include <string>
#include <unordered_set>
#include <cassert>
#include "dag_task.h"


// An abstract representation for schedules of subtasks.
// A subtask may be scheduled into multiple fragments.
struct Fragment {
	// ID of this specific fragment.
	uint_t frag_id;

	// ID of the corresponding subtask. A subtask can have more than 1 fragments.
	uint_t sub_id;
	
	// Time when the fragment is executed.
	uint_t start_time;

	// Execution time of this fragment.
	uint_t exe_time;

	// The time by which this fragment must finish.
	uint_t deadline;
	
	// The total work of subdag starting from this fragment.
	// The subdag includes all fragments and nodes up to some sink nodes.
	// We use this field to compute the corresponding density of the
	// subdag and use that density to decide which fragment is scheduled
	// next.
	uint_t subdag_work;
	
	// Processor ID on which this fragment executes. IDs start from 0.
	uint_t proc_id;

	// Split this fragment into 2 fragments: one is being scheduled,
	// the other is the remaining part waiting to be scheduled.
	// Input: start time and execution time of the scheduled part, and
	// ID of the processor on which it is scheduled.
	// Output: the remaining part of the fragments.
	// Note that this method assumes @exe_time is less than the current
	// execution time of the split fragment.
	Fragment split(uint_t start, uint_t work, uint_t proc_id) {
		// Make sure the assumption holds.
		assert(this->exe_time > work);

		// The remaining fragment.
		Fragment remain = *this;
		remain.frag_id++;
		remain.exe_time -= work;
		remain.subdag_work -= work;
		
		// Update the first part.
		this->start_time = start;
		this->exe_time = work;
		this->proc_id = proc_id;
		return remain;
	}

	// Set the start time and processor ID for the last fragment.
	void update_final_fragment(uint_t start, uint_t proc_id) {
		this->start_time = start;
		this->proc_id = proc_id;
	}
};


// Subtask can have multiple fragments. 
struct Subtask {
	uint_t sub_id;
	list<Fragment> frags;
};

// The length of a longest path starting from a fragment of a subtask.
struct LongestPath {
	uint_t sub_id;
	uint_t length;
};

// Compare the longest-path lengths of starting from the corresponding
// fragments. The one with greatest length is placed on top.
struct GreaterLength {
	bool operator()(const LongestPath& a, const LongestPath& b) const {
		return a.length < b.length;
	}
};

// Used to compare the works of sub-dags starting from the
// corresponding fragments. The one has greatest work is on top.
struct GreaterDensity {
	bool operator()(const Fragment& a, const Fragment& b) const {
		return a.subdag_work < b.subdag_work;
	}
};

class ProposedScheduler {
public:
	ProposedScheduler(const DagTask& t);

	const vector<Subtask>& get_schedule() const {
		return subtasks_;
	}

	uint_t get_num_procs() const {
		return num_procs_;
	}

	void write_task_and_schedule(string, string) const;
	void write_schedule(string) const;
	void write_dagtask(string) const;
	
private:
	void initialize();
	void schedule();
	bool schedule_core();
	void add_ready_fragments(uint_t, priority_queue<Fragment, vector<Fragment>, GreaterDensity>&, unordered_set<uint_t>&);
	bool schedule_core_old();
	void add_ready_fragments_old(uint_t, priority_queue<Fragment, vector<Fragment>, GreaterDensity>&);
	double density(uint_t, uint_t);
	
private:
	DagTask task_;
	
	// This list is updated while the scheduling procedure runs,
	// so we have to get another copy of it when rerun the procedure.
	list_t depend_list_;
	
	// Number of processors for this task.
	uint_t num_procs_;
	
	// Update the current time of the schedule while constructing it.
	uint_t curr_time_;
	uint_t deadline_;

	// Keep the remaining work of the DAG while constructing its schedule.
	uint_t work_left_;
	
	// Update the schedule for the fragments while the algorithm runs.
	vector<Subtask> subtasks_;
};

#endif
