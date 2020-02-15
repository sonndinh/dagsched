#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "naive_greedy_scheduler.h"

NaiveGreedyScheduler::NaiveGreedyScheduler(const DagTask& t) : task_{t} {
	uint_t work = task_.get_work();
	uint_t span = task_.get_span();
	uint_t deadline = task_.get_deadline();

	// Be careful with the edge cases when deadline equals span OR
	// work equals span.
	if (work == span) {
		// This task can only execute sequentially. So it is only
		// allocated 1 processor.
		num_procs_ = 1;
	} else if (deadline == span) {
		// This is meant to say the schedule needs the number of
		// processors equal to its maximum parallelism.
		num_procs_ = numeric_limits<unsigned int>::max();
	} else {
		// Normal cases.
		num_procs_ = ceil((double)(work - span)/(deadline - span));
	}
	if (num_procs_ < 1) {
		cout << "Work: " << work << ". Span: " << span << ". Deadline: " << deadline << ". #procs: " << num_procs_ << endl;
	}
	
	depend_list_ = task_.get_depend_list();
	if (num_procs_ < numeric_limits<unsigned int>::max()) {
		procs_.resize(num_procs_);
		for (int i = 0; i < num_procs_; ++i) {
			procs_[i].id = i;
		}
	}
	
	// Schedule the DAG task.
	schedule();
}

NaiveGreedyScheduler::NaiveGreedyScheduler() : num_procs_{0} {
}

void NaiveGreedyScheduler::unrestricted_schedule(const DagTask& t) {
	task_ = t;
	depend_list_ = task_.get_depend_list();
	build_trivial_schedule();
}

// Trivial schedule for the task. Called when deadline is equal to span.
// Basically it is similar to the schedule() method except that we add
// new processors when needed while the method runs. 
void NaiveGreedyScheduler::build_trivial_schedule() {
	priority_queue<Processor, vector<Processor>, AvailTimeCompare> free_procs;
	// Insert one processor.
	free_procs.push(Processor());
	const vector<uint_t>& sources = task_.get_sources();
	priority_queue<ReadySubtask, vector<ReadySubtask>, ReadyTimeCompare> ready_subtasks;
	// Keep the start times for all subtasks. Each equals to the largest finish time of its parents.
	vector<uint_t> start_times(task_.get_size(), 0);
	
	for (uint_t id : sources) {
		ready_subtasks.push({id, 0});
	}

	uint_t count_procs = 1;
	while (!ready_subtasks.empty()) {
		ReadySubtask next_subtask = ready_subtasks.top();
		uint_t sub_id = next_subtask.id;
		ready_subtasks.pop();
		uint_t subtask_work = task_.get_node_lengths()[sub_id];
		Processor proc = free_procs.top();
		
		if (proc.avail_time > next_subtask.ready_time) {
			// Add a new processor and schedule the next subtask on it.
			Processor new_proc;
			new_proc.id = count_procs++;
			new_proc.avail_time = next_subtask.ready_time + subtask_work;
			new_proc.subtasks.push_back({sub_id, next_subtask.ready_time});
			free_procs.push(new_proc);
		} else {
			// Run this subtask on the existing processor.
			free_procs.pop();
			proc.avail_time = next_subtask.ready_time + subtask_work;
			proc.subtasks.push_back({sub_id, next_subtask.ready_time});
			free_procs.push(proc);
		}

		// Add new ready subtasks.
		add_ready_subtasks(ready_subtasks, start_times, sub_id, next_subtask.ready_time + subtask_work);
	}
	
	// Obtain the schedule of the processors.
	procs_.resize(free_procs.size());
	while (!free_procs.empty()) {
		Processor proc = free_procs.top();
		free_procs.pop();
		procs_[proc.id] = proc;
	}
}


// Return the finish time of the whole DAG task.
uint_t NaiveGreedyScheduler::schedule() {
	if (num_procs_ == numeric_limits<unsigned int>::max()) {
		// This happens when the deadline equals to the span of the task.
		// Trivial schedule constructed by simply taking the DAG structure.
		// The task finishes exactly at its deadline.
		build_trivial_schedule();
		return task_.get_deadline();
	}
		
	// Top entry is processor with earliest available time.
	priority_queue<Processor, vector<Processor>, AvailTimeCompare> free_procs(procs_.begin(), procs_.end());
	const vector<uint_t>& sources = task_.get_sources();
	priority_queue<ReadySubtask, vector<ReadySubtask>, ReadyTimeCompare> ready_subtasks;
	// Keep the start times for all subtasks. Each equals to the largest finish time of its parents.
	vector<uint_t> start_times(task_.get_size(), 0);
	// Keep the latest finish time of the subtasks.
	uint_t finish_time = 0;

	for (uint_t id : sources) {
		ready_subtasks.push({id, 0});
	}

	// In each iteration, schedule the subtask with earliest ready time on
	// the processor with earliest available time.
	// Add new ready subtasks to the priority queue.
	while (!ready_subtasks.empty()) {
		ReadySubtask next_subtask = ready_subtasks.top();
		uint_t sub_id = next_subtask.id;
		ready_subtasks.pop();
		// Debugging -- start
		//if (free_procs.empty()) {
		//	cout << "ERROR: Empty queue for available processors !!" << endl;
		//}
		// Debugging -- end
		Processor proc = free_procs.top();
		free_procs.pop();
		uint_t subtask_work = task_.get_node_lengths()[sub_id];

		// Update the new available time for the chosen processor.
		uint_t start_time = max(next_subtask.ready_time, proc.avail_time);
		proc.avail_time = start_time + subtask_work;
		proc.subtasks.push_back({sub_id, start_time});
		finish_time = max(finish_time, proc.avail_time);

		// Add new ready subtasks.
		add_ready_subtasks(ready_subtasks, start_times, sub_id, start_time + subtask_work);

		// Push the processor back to the heap.
		free_procs.push(proc);
	}

	// Obtain the schedule of the processors.
	while (!free_procs.empty()) {
		Processor proc = free_procs.top();
		free_procs.pop();
		procs_[proc.id] = proc;
	}

	//	cout << "Finish time: " << finish_time << endl;
	return finish_time;
}


// Add new ready subtasks after subtask @id finishes.
//void NaiveGreedyScheduler::add_ready_subtasks(queue<uint_t>& ready_subtasks, uint_t id) {
void NaiveGreedyScheduler::add_ready_subtasks(priority_queue<ReadySubtask, vector<ReadySubtask>, ReadyTimeCompare>& ready_subtasks,
											  vector<uint_t>& start_times, uint_t id, uint_t finished_time) {
	const list<uint_t>& adj = task_.get_adj_list()[id];
	for (list<uint_t>::const_iterator it = adj.cbegin(); it != adj.cend(); ++it) {
		uint_t next_sub_id = *it;
		// Update the earliest start time possible for each adjacent subtask.
		start_times[next_sub_id] = max(start_times[next_sub_id], finished_time);
		
		list<uint_t>& parents = depend_list_[next_sub_id];
		list<uint_t>::const_iterator pos = find(parents.cbegin(), parents.cend(), id);
		parents.erase(pos);

		// All parents have finished. This subtask is ready now.
		if (parents.empty()) {
			ready_subtasks.push({next_sub_id, start_times[next_sub_id]});
		}
	}
}

// Write the schedule to the given file.
void NaiveGreedyScheduler::write_schedule(string file_path) const {
	fstream fs(file_path, ios_base::out);
	if (!fs.is_open()) {
		cout << "ERROR: Could not open file " << file_path << " !!" << endl;
		return ;
	}

	const vector<uint_t>& wcets = task_.get_node_lengths();
	
	// First line is the number of processors assigned to this task.
	fs << procs_.size() << endl;

	// Each consequent line contains a list of subtasks scheduled on
	// the corresponding processor: <processor id> followed by the list 
	// of subtasks <subtask id> <start time> <execution time>.
	for (const Processor& proc : procs_) {
		fs << proc.id << " ";
		for (const FinishedSubtask& sub : proc.subtasks) {
			fs << sub.id << " " << sub.start_time << " " << wcets[sub.id] << " ";
		}
		fs << endl;
	}
	
	fs.close();
}
