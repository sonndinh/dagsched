#include <iostream>
#include <cassert>
#include <limits>
#include <cmath>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "proposed_scheduler.h"

ProposedScheduler::ProposedScheduler(const DagTask& t) : task_{t}, num_procs_{1}, curr_time_{0} {
	// Initialize the fragments data.
	initialize();

	// Schedule the task.
	schedule();
}

// Initialize the data at the beginning of the scheduling procedure.
void ProposedScheduler::initialize() {
	// Reset the data.
	curr_time_ = 0;
	depend_list_ = task_.get_depend_list();
	deadline_ = task_.get_deadline();
	work_left_ = task_.get_work();	
	
	uint_t size = task_.get_size();
	const vector<uint_t>& sources = task_.get_sources();
	const vector<uint_t>& sinks = task_.get_sinks();
	subtasks_.resize(size);

	// Init original fragments corresponding to the subtasks.
	// At begin, each fragment is considered as a whole subtask.
	for (int i = 0; i < size; ++i) {
		// Clear the old list of fragments.
		subtasks_[i].frags.clear();
		subtasks_[i].sub_id = i;
		Fragment frag;
		frag.frag_id = 0;
		frag.sub_id = i;

		// Random initial value for @start_time.
		frag.start_time = numeric_limits<uint_t>::max();
		
		frag.exe_time = task_.get_node_lengths()[i];

		// Initialize the deadline to the DAG task's deadline.
		// Update it as the scheduling algorithm runs.
		frag.deadline = deadline_;
		
		// Set the fragment's subdag work.
		frag.subdag_work = task_.get_subdag_works()[i];
		
		// Random initial value for @proc_id.d
		frag.proc_id = numeric_limits<uint_t>::max();
		
		subtasks_[i].frags.push_back(frag);
	}
}

// Constructor calls this method.
void ProposedScheduler::schedule() {
	while (schedule_core() == false) {
		// Reset data everytime we run the scheduling procedure.
		initialize();
	}
}

double ProposedScheduler::density(uint_t work, uint_t deadline) {
	return (double)work/deadline;
}

// Core schedule method. Return true if it schedules successfully.
// Return false if it needs to run again with a new number of processors.
bool ProposedScheduler::schedule_core() {
	const vector<uint_t> &sources = task_.get_sources();
	// A ready queue that contains all ready fragments.
	priority_queue<Fragment, vector<Fragment>, GreaterDensity> density_q;
	// Lengths of longest paths starting from every subtasks.
	vector<uint_t> path_lengths = task_.get_longest_paths();
	// This set maintains the IDs of ready subtasks while the while-loop runs.
	// It contains IDs for the same subtasks as in density_q.
	// We must keep this set and density_q consistent with each other during the loop.
	unordered_set<uint_t> ready_sub_ids;

	// Initialize the two priority queues.
	for (uint_t id : sources) {
		density_q.push(subtasks_[id].frags.front());
		ready_sub_ids.insert(id);
	}

	while (!density_q.empty()) {
		// There are 2 cases when we need to rerun the scheduling procedure.
		// 1. When the lower-bound on the number of processors required to schedule
		// the remaining work is greater than the current number of processors.
		// 2. When there are more paths that need to be executed immediately, i.e.,
		// paths with lengths equal to the time left until deadline, than the 
		// current number of processors. 
		uint_t estimated_num_procs = ceil(density(work_left_, deadline_ - curr_time_));
		if (estimated_num_procs > num_procs_) {
			// Only increase the number of processors by 1 in this version.
			// This more conservative increment is likely better.
			++num_procs_;
			return false;
		}
		
		// Number of paths whose lengths equal to the time left until deadline.
		uint_t count_full_paths = 0;
		uint_t time_left = deadline_ - curr_time_;
		for (uint_t id : ready_sub_ids) {
			if (path_lengths[id] > time_left) {
				cout << "======= ERROR: A path of subtask " << id << " will miss deadline!!" << endl;
			}
			if (path_lengths[id] == time_left) {
				++count_full_paths;
			}
		}
		if (count_full_paths > num_procs_) {
			// Cannot schedule successfully with the current number of processors.
			// So we increase the no. of processors by 1 and rerun the procedure.
			++num_procs_;
			return false;
		}

		// At this point, we can continue with the current number of processors.
		// Pick at most num_procs_ ready fragments to execute. First, we pick all ready fragments
		// with longest path length equal to time_left, i.e., those paths that will miss deadline
		// if not executed immediately. After that, we pick among the remaining ready fragments
		// the ones with greatest densities up to total of num_procs_ fragments to run in this iteration.
		vector<Fragment> frags_being_scheduled;
		// IDs of subtasks corresponding to the longest paths that need to run immediately.
		unordered_set<uint_t> frags_from_longest_paths;
		uint_t count = 0; // The number of chosen fragments.
		for (uint_t id : ready_sub_ids) {
			if (path_lengths[id] == time_left) {
				frags_from_longest_paths.insert(id);
				++count;
			}
		}

		// Add the ready fragments corresponding to those longest paths.
		uint_t i = 0;
		vector<Fragment> temp;
		while (!density_q.empty() && i < count) {
			uint_t id = density_q.top().sub_id;
			if (frags_from_longest_paths.find(id) != frags_from_longest_paths.cend()) {
				frags_being_scheduled.push_back(density_q.top());
				density_q.pop();
				ready_sub_ids.erase(id); // Remove from the ready set too.
				++i;
			} else {
				// Temporarily pop this fragment out so we can examine the next fragments.
				temp.push_back(density_q.top());
				density_q.pop();
			}
		}
		// Push back those fragments from temp to the ready queue.
		for (Fragment &frag : temp) {
			density_q.push(frag);
		}
				
		// Pick among the remaining fragments the ones with highest densities (i.e., greatest subgraph work).
		// The smallest subgraph work among the scheduled fragments.
		uint_t smallest_subdag_work = numeric_limits<unsigned int>::max();
		while (!density_q.empty() && count < num_procs_) {
			frags_being_scheduled.push_back(density_q.top());
			if (density_q.top().subdag_work < smallest_subdag_work) {
				smallest_subdag_work = density_q.top().subdag_work;
			}
			ready_sub_ids.erase(density_q.top().sub_id); // Remove from the ready set too.
			density_q.pop();
			++count;
		}

		// The largest subdag work from unscheduled fragments in this iteration.
		uint_t unsched_max_work = 0;
		if (!density_q.empty()) {
			unsched_max_work = density_q.top().subdag_work;
		}

		// Compute the execution length for those fragments.
		uint_t min_exec_time = numeric_limits<unsigned int>::max();
		for (auto &frag : frags_being_scheduled) {
			min_exec_time = min(min_exec_time, frag.exe_time);
		}

		// Find the length of a ready fragment with longest path among 
		// those which are not scheduled in this iteration.
		uint_t unsched_max_length = 0;
		for (uint_t id : ready_sub_ids) {
			bool scheduled = false;
			for (auto &frag : frags_being_scheduled) {
				if (id == frag.sub_id) {
					scheduled = true;
					break;
				}
			}
			if (scheduled == false) {
				unsched_max_length = max(unsched_max_length, path_lengths[id]);
			}
		}
		if (time_left < unsched_max_length) {
			// Should NEVER happen!
			cout << "======= ERROR: Unscheduled path length: " << unsched_max_length
				 << " < Time left: " << time_left << " !!" << endl;
		}
		min_exec_time = min(min_exec_time, time_left - unsched_max_length);

		// We compute the difference (denoted by Delta) between the smallest total
		// work among the fragments chosen to be scheduled in this iteration and
		// the largest total work among the fragments not chosen to be scheduled.
		// The execution time for the chosen fragments in this iteration will be at most: Delta + 1.
		// The reason for using Delta is that we only want to execute the chosen fragments until
		// the time when one of the unscheduled fragments has remaining work larger than
		// the remaining work of at least one of the chosen fragments. At that point, we want to
		// make a scheduling decision again. 1 unit is added to Delta to break ties in case there are
		// multiple fragments with the same total remaining work.
		if (smallest_subdag_work < numeric_limits<unsigned int>::max() && unsched_max_work > 0) {
			min_exec_time = min(min_exec_time, smallest_subdag_work - unsched_max_work + 1);
		}
		
		// Debugging -- start
#ifdef _DEBUG_
		cout << "======== Number of processors: " << num_procs_ << endl;
		cout << "t=" << curr_time_ << ": Run ";
		for (int i = 0; i < frags_being_scheduled.size(); ++i) {
			uint_t sub_id = frags_being_scheduled[i].sub_id;
			uint_t frag_id = frags_being_scheduled[i].frag_id;
			cout << "(subtask " << sub_id << ",frag " << frag_id << ") ";
		}
		cout << " for " << min_exec_time << " units." << endl;
#endif
		// Debugging -- end

		// If the previous fragment of a chosen fragment was scheduled in the immediately
		// preceding iteration, schedule the chosen fragment on the same processor.
		// Set of processors on which fragments can be schedule. 
		unordered_set<uint_t> free_procs;
		for (uint_t i = 0; i < num_procs_; ++i) {
			free_procs.insert(i);
		}
		// Map subtasks from previous iteration to their corresponding processors. 
		unordered_map<uint_t, uint_t> sub_to_proc;
		for (int i = 0; i < frags_being_scheduled.size(); ++i) {
			uint_t sub_id = frags_being_scheduled[i].sub_id;
			if (subtasks_[sub_id].frags.size() < 2)
				continue;
			
			auto previous = --subtasks_[sub_id].frags.end();
			--previous;
			uint_t finish_time = previous->start_time + previous->exe_time;
			if (finish_time == curr_time_) {
				sub_to_proc[sub_id] = previous->proc_id;
				free_procs.erase(previous->proc_id);
			}
		}
		
		// Split the fragments, enable new ready subtasks, update tracking data.
		// Schedule those fragments and split them if necessary.
		for (int i = 0; i < frags_being_scheduled.size(); ++i) {
			uint_t sub_id = frags_being_scheduled[i].sub_id;
			uint_t frag_id = frags_being_scheduled[i].frag_id;
			Subtask& subtask = subtasks_[sub_id];
			
			// If the scheduled fragment is the first fragment of its subtask,
			// update the deadline of the last fragments of its parent subtasks.
			// Otherwise, update the deadline of its immediate prior fragment from the same subtask.
			if (subtask.frags.size() == 1) {
				const list<uint_t>& parents = task_.get_depend_list()[sub_id];
				for (uint_t parent_id : parents) {
					uint_t curr_deadline = subtasks_[parent_id].frags.back().deadline;
					subtasks_[parent_id].frags.back().deadline = min(curr_deadline, curr_time_);
				}
			} else {
				list<Fragment>::iterator prev = --subtask.frags.end();
				--prev;
				prev->deadline = curr_time_;
			}

			// Get appropriate processor for this fragment.
			uint_t proc_id;
			if (sub_to_proc.find(sub_id) != sub_to_proc.end()) {
				proc_id = sub_to_proc[sub_id];
			} else { // Assign to any other processor.
				proc_id = *(free_procs.begin());
				free_procs.erase(free_procs.begin());
			}
			
			if (subtask.frags.back().exe_time > min_exec_time) {
				// Split this fragment.
				Fragment remain = subtask.frags.back().split(curr_time_, min_exec_time, proc_id);
				density_q.push(remain);
				ready_sub_ids.insert(remain.sub_id); // Add the subtask ID back to the set.
				subtask.frags.push_back(remain);
			} else {
				// This fragment is scheduled as a whole and
				// its corresponding subtask is finished.
				subtask.frags.back().update_final_fragment(curr_time_, proc_id);
				
				// Enable the depending subtasks and add them to ready queue and set.
				add_ready_fragments(sub_id, density_q, ready_sub_ids);
			}
		}		

		// Update path lengths for the paths that are scheduled.
		for (auto &frag : frags_being_scheduled) {
			uint_t id = frag.sub_id;
			path_lengths[id] -= min_exec_time;
		}
		
		// Advance the current time.
		curr_time_ += min_exec_time;
		
		// Update the amount of work left.
		work_left_ -= frags_being_scheduled.size() * min_exec_time;
	}
	
	return true;
}

// Add new enabled subtasks to the ready queue and the set of IDs for ready subtasks.
void ProposedScheduler::add_ready_fragments(uint_t finished_sub_id, priority_queue<Fragment, vector<Fragment>, GreaterDensity>& ready_frags,
											unordered_set<uint_t>& ready_ids) {
	const list<uint_t>& adj = task_.get_adj_list()[finished_sub_id];
	for (list<uint_t>::const_iterator it = adj.cbegin(); it != adj.cend(); ++it) {
		uint_t next_sub_id = *it;
		list<uint_t>& parents = depend_list_[next_sub_id];
		list<uint_t>::const_iterator pos = find(parents.cbegin(), parents.cend(), finished_sub_id);
		// Note that the dependency list is modified here.
		parents.erase(pos);
		
		if (parents.empty()) {
			// This subtask is ready to execute now.
			ready_frags.push(subtasks_[next_sub_id].frags.back());
			ready_ids.insert(next_sub_id);
		}
	}
}

// DEPRECATED version of schedule_core() method.
// There are tasks for which the schedules computed by this method
// miss their deadlines. This is because this method does NOT consider
// the lengths of the longest paths starting from the ready fragments.
bool ProposedScheduler::schedule_core_old() {
	vector<uint_t> crit_path = task_.get_crit_path();
	uint_t slack = deadline_ - task_.get_span();
	const vector<uint_t>& sources = task_.get_sources();
	uint_t span_left = task_.get_span(); // Remaining length of critical-path.
	priority_queue<Fragment, vector<Fragment>, GreaterDensity> ready_frags;

	// Add source subtasks.
	for (uint_t id : sources) {
		ready_frags.push(subtasks_[id].frags.front());
	}

	while (!ready_frags.empty()) {
		uint_t estimated_num_procs = ceil(density(work_left_, deadline_ - curr_time_));
		if (estimated_num_procs > num_procs_) {
			// For simplicity, whenever the estimated number of processors required	increases,
			// we restart the scheduling method using the new estimated number of processors.
			num_procs_ = estimated_num_procs;
			return false;
		}
		
		// Pop out num_procs ready fragments.
		vector<Fragment> frags_being_scheduled;
		int i = 0;
		bool crit_frag_left = false;
		Fragment crit_frag;
		while (!ready_frags.empty() && i < num_procs_) {
			Fragment frag = ready_frags.top();
			ready_frags.pop();
			if (find(crit_path.begin(), crit_path.end(), frag.sub_id) == crit_path.end()) {
				if (crit_frag_left) {
					double crit_density = density(crit_frag.subdag_work, deadline_ - curr_time_);
					// Gotta schedule the critical fragment if either there is no more slack OR
					// the critical fragment's density is greater than the next fragment's and 1.0.
					if (slack == 0 || (frag.subdag_work < crit_frag.subdag_work && crit_density >= 1.0)) {
						frags_being_scheduled.push_back(crit_frag);
						crit_frag_left = false;
						++i;
					}
				}
				if (i < num_procs_) {
					frags_being_scheduled.push_back(frag);
					++i;
				}
			} else {
				// The next fragment with largest density is a critical one.
				// Postpone the decision of whether to schedule this fragment to the next iteration.
				crit_frag_left = true;
				crit_frag = frag;
			}
		}
		if (crit_frag_left) {
			if (i < num_procs_) {
				// There is a processor for scheduling the critical fragment.
				frags_being_scheduled.push_back(crit_frag);
			} else {
				// Otherwise, push it back to the queue for considering in the next round.
				ready_frags.push(crit_frag);
			}
		}
		
		// Find the minimum execution time among those fragments.
		uint_t min_exec_time = numeric_limits<unsigned int>::max();
		bool has_critical_fragment = false;
		for (int i = 0; i < frags_being_scheduled.size(); ++i) {
			min_exec_time = min(min_exec_time, frags_being_scheduled[i].exe_time);
			uint_t sub_id = frags_being_scheduled[i].sub_id;
			if (find(crit_path.begin(), crit_path.end(), sub_id) != crit_path.end()) {
				has_critical_fragment = true;
			}
		}
		
 		if (span_left > 0 && slack == 0 && !has_critical_fragment) {
			// The critical-path has NOT been scheduled completely yet and there is no more slack left,
			// but the scheduler did not choose to schedule a critical fragment in this step.
			// This will lead to the critical-path to run pass the deadline!
			cout << "FATAL: Current time: " << curr_time_ << ". The critical-path will run pass deadline !!" << endl;
		}
		
		if (slack > 0 && !has_critical_fragment) {
			// Update the slack when no critical fragment is scheduled in this step.
			if (slack > min_exec_time) {
				slack -= min_exec_time;
			} else {
				min_exec_time = slack;
				slack = 0;
			}
		}

		// Schedule those fragments and split them if necessary.
		for (int i = 0; i < frags_being_scheduled.size(); ++i) {
			uint_t sub_id = frags_being_scheduled[i].sub_id;
			uint_t frag_id = frags_being_scheduled[i].frag_id;
			Subtask& subtask = subtasks_[sub_id];

			// If the scheduled fragment is the first fragment of its subtask,
			// update the deadline of the last fragments of its parent subtasks.
			// Otherwise, update the deadline of its immediate prior fragment from the same subtask.
			if (subtask.frags.size() == 1) {
				const list<uint_t>& parents = task_.get_depend_list()[sub_id];
				for (uint_t parent_id : parents) {
					uint_t curr_deadline = subtasks_[parent_id].frags.back().deadline;
					subtasks_[parent_id].frags.back().deadline = min(curr_deadline, curr_time_);
				}
			} else {
				list<Fragment>::iterator prev = --subtask.frags.end();
				--prev;
				prev->deadline = curr_time_;
			}

			
			if (subtask.frags.back().exe_time > min_exec_time) {
				// Split this fragment.
				Fragment remain = subtask.frags.back().split(curr_time_, min_exec_time, i);
				ready_frags.push(remain);
				subtask.frags.push_back(remain);
			} else {
				// This fragment is scheduled as a whole and
				// its corresponding subtask is finished.
				subtask.frags.back().update_final_fragment(curr_time_, i);
				
				// Enable the depending subtasks and add to ready queue.
				add_ready_fragments_old(sub_id, ready_frags);
			}
		}

		// Advance the current time.
		curr_time_ += min_exec_time;
		
		// Update the amount of work left.
		work_left_ -= frags_being_scheduled.size() * min_exec_time;

		if (has_critical_fragment) {
			// A critical fragment was scheduled. So we need to
			// update the remaining length of the critical-path.
			span_left -= min_exec_time;
		}
	}
	return true;
}


void ProposedScheduler::add_ready_fragments_old(uint_t finished_sub_id, priority_queue<Fragment, vector<Fragment>, GreaterDensity>& ready_frags) {
	const list<uint_t>& adj = task_.get_adj_list()[finished_sub_id];
	for (list<uint_t>::const_iterator it = adj.cbegin(); it != adj.cend(); ++it) {
		uint_t next_sub_id = *it;
		list<uint_t>& parents = depend_list_[next_sub_id];
		list<uint_t>::const_iterator pos = find(parents.cbegin(), parents.cend(), finished_sub_id);
		// Note that the dependency list is modified here.
		parents.erase(pos);
		
		if (parents.empty()) {
			// This subtask is ready to execute now.
			ready_frags.push(subtasks_[next_sub_id].frags.back());
		}
	}
}

// Write both task structure and its resulting schedule.
void ProposedScheduler::write_task_and_schedule(string task_file, string schedule_file) const {
	write_dagtask(task_file);
	write_schedule(schedule_file);
}


// Write the obtained schedule to a file with the given path.
void ProposedScheduler::write_schedule(string file_path) const {
	fstream fs(file_path, ios_base::out);
	if (!fs.is_open()) {
		cout << "ERROR: Could not open file " << file_path << " !!" << endl;
		return;
	}

	// First line is the number of subtasks.
	fs << subtasks_.size() << endl;

	// Sep 03, 2019: Add a second line for the number of processors and
	// and the finish time for the DAG task.
	fs << num_procs_ << " " << curr_time_ << endl;

	// Each subsequent line contains the list of fragments for
	// the corresponding subtask.
	// Each contains: <subtask id> <number of fragments>, then for each
	// fragment: <processor id> <start time> <execution time> <deadline>
	for (const Subtask& s : subtasks_) {
		fs << s.sub_id << " " << s.frags.size() << " ";
		for (const Fragment& frag : s.frags) {
			fs << frag.proc_id << " " << frag.start_time << " "
				 << frag.exe_time << " " << frag.deadline << " ";
		}
		fs << endl;
	}
	
	fs.close();
}


// Write DAG task's adjacency list to the given file.
void ProposedScheduler::write_dagtask(string file_path) const {
	fstream fs(file_path, ios_base::out);
	if (!fs.is_open()) {
		cout << "ERROR: Could not open file " << file_path << " !!" << endl;
		return;
	}

	const list_t& adj_list = task_.get_adj_list();

	// First line is the number of subtasks.
	fs << adj_list.size() << endl;

	// Second line is the WCETs of the subtasks.
	auto wcets = task_.get_node_lengths();
	for (uint_t node_len : wcets) {
		fs << node_len << " ";
	}
	fs << endl;

	// Third line contains the work, span, deadline, period.
	// Note: Sep 03, 2019: Added work to this line.
	fs << task_.get_work() << " " << task_.get_span() << " "
	   << task_.get_deadline() << " " << task_.get_period() << endl;

	// Forth line contains IDs of subtasks on the critical-path.
	const vector<uint_t>& crit_path = task_.get_crit_path();
	for (uint_t sub_id : crit_path) {
		fs << sub_id << " ";
	}
	fs << endl;

	// Each subsequent line contains a list of adjacent subtasks of the
	// corresponding subtask: <subtask id> followed by ids of the adjacent subtasks.
	uint_t sub_id = 0;
	for (const list<uint_t>& adj : adj_list) {
		fs << sub_id << " ";
		++sub_id;
		for (uint_t neighbor_id : adj) {
			fs << neighbor_id << " ";
		}
		fs << endl;
	}
	
	fs.close();
}
