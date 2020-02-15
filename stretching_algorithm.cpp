#include <limits>
#include <algorithm>
#include <queue>
#include "stretching_algorithm.h"
#include "naive_greedy_scheduler.h"

StretchingAlgo::StretchingAlgo(const Taskset& t, SchedulingAlgo algo) : tset_{t}, algo_{algo} {}

// Build Multi-Threaded Segment representation for the input task. 
void StretchingAlgo::build_mts(const DagTask& task, const NaiveGreedyScheduler& sched) {
	// List of segments for this particular task.
	mts_t mts;
	vector<Processor> procs_sched = sched.get_schedule();
	vector<int> indexes(procs_sched.size(), 0);
	uint_t curr_time = 0;
	uint_t deadline = task.get_deadline();

	// Copy the schedule to vector for easier access.
	vector<vector<FinishedSubtask>> procs;
	for (const Processor& p : procs_sched) {
		procs.emplace_back(vector<FinishedSubtask>(p.subtasks.begin(), p.subtasks.end()));
	}

	while (true) {
		uint_t segment_len = numeric_limits<unsigned int>::max();
		Segment seg;
		// Compute the length of the next segment.
		for (uint_t i = 0; i < procs.size(); ++i) {
			int idx = indexes[i];
			if (idx == -1)
				continue;
			
			uint_t node_id = procs[i][idx].id;
			uint_t node_start = procs[i][idx].start_time;
			uint_t node_len = task.get_node_lengths()[node_id];
			uint_t node_end = node_start + node_len;
			
			if (node_end > curr_time && node_start <= curr_time) {
				segment_len = min(segment_len, node_end - curr_time);
				seg.subtasks.push_back(node_id);
			}
		}
		if (segment_len == numeric_limits<unsigned int>::max())
			break;

		// Add a new segment and update bookkeeping variables.
		seg.start = curr_time;
		seg.len = segment_len;
		mts.push_back(seg);
		
		curr_time += segment_len;
		for (uint_t i = 0; i < procs.size(); ++i) {
			int idx = indexes[i];
			if (idx == -1)
				continue;

			uint_t node_id = procs[i][idx].id;
			uint_t node_start = procs[i][idx].start_time;
			uint_t node_len = task.get_node_lengths()[node_id];
			uint_t node_end = node_start + node_len;
			if (node_end == curr_time) {
				if (idx == procs[i].size() - 1)
					indexes[i] = -1;
				else
					indexes[i] += 1;
			}
		}
	}

	// Add the MTS of this task.
	mtss_.push_back(mts);
}

// Stretch this heavy task.
void StretchingAlgo::stretch(uint_t id, mts_t& mts) {
	uint_t work = tset_.get_heavy_tasks()[id].get_work();
	uint_t span = tset_.get_heavy_tasks()[id].get_span();
	uint_t deadline = tset_.get_heavy_tasks()[id].get_deadline();
	uint_t period = tset_.get_heavy_tasks()[id].get_period();
	
	// Distribution factor.
	double f = (double)(deadline - span)/(work - span);

	// Set of the threads remained after stretching this DAG task. 
	vector<SequentialThread> tmp;
	
	// Go through each segment and stretch it.
	for (uint_t i = 0; i < mts.size(); ++i) {
		Segment& seg = mts[i];
		// Number of additional threads being added to the master
		// thread, not including the critical thread. 
		double f_j = f * (seg.subtasks.size() - 1);
		// Number of full threads added to the master thread.
		uint_t q_j = (uint_t)floor(f_j) + 1;
		
		// Update the start time for the segment if it is not the first one.
		// The first segment still has the start time of 0.
		if (i > 0) {
			// The length of the (i-1)th segment was updated after stretching.
			seg.start = mts[i-1].start + mts[i-1].len;
		}

		if (floor(f_j) < f_j) {
			// The first part of a split thread.
			double frac_work = (1 + floor(f_j) - f_j) * seg.len;
			double frac_deadline = (1 + floor(f_j)) * seg.len;
			// Add the first part of the split thread to set of sequential tasks.
			tmp.push_back({seg.start, frac_work, frac_deadline, period});
		}
		
		// Add the other threads.
		for (uint_t k = q_j + 2; k <= seg.subtasks.size(); ++k) {
			double my_work = seg.len;
			double my_deadline = (1 + f_j) * seg.len;
			tmp.push_back({seg.start, my_work, my_deadline, period});
		}

		// Update the length of this segment.
		seg.len *= (1 + f_j);
	}

	// Add the sequential threads for this DAG task.
	threads_.push_back(tmp);
}

bool StretchingAlgo::is_schedulable() {
	// Compute the schedule for each DAG task on an unrestricted number of
	// processors, and convert to the Multi-Threaded Segment (MTS) representation.
	for (const DagTask& t : tset_.get_heavy_tasks()) {
		NaiveGreedyScheduler naive_sched;
		naive_sched.unrestricted_schedule(t);
		build_mts(t, naive_sched);
	}
	// Testing...
	//	print_mts();

	// Stretch the heavy DAG tasks.
	for (uint_t i = 0; i < mtss_.size(); ++i) {
		stretch(i, mtss_[i]);
	}
	// Testing...
	//	print_stretched_threads();

	// Number of dedicated cores for master threads. 
	uint_t required_procs = tset_.get_heavy_tasks().size();
	//cout << "Number of dedicated cores for heavy tasks: " << required_procs << endl;
	
	// Add the threads of light DAG tasks.
	for (uint_t i = 0; i < tset_.get_light_tasks().size(); ++i) {
		double work = tset_.get_light_tasks()[i].get_work();
		double deadline = tset_.get_light_tasks()[i].get_deadline();
		uint_t period = tset_.get_light_tasks()[i].get_period();
		if (work == deadline) {
			++required_procs;
		} else {
			vector<SequentialThread> tmp;
			tmp.push_back({0, work, deadline, period});
			threads_.push_back(tmp);
		}
	}
	//cout << "Total number of dedicated cores: " << required_procs << endl;
	
	if (tset_.get_no_procs() < required_procs) {
		//cout << "Total no. cores: " << tset_.get_no_procs() << " < required cores !!" << endl;
		return false;
	}

	uint_t remain_cores = tset_.get_no_procs() - required_procs;
	// Each master thread gets one dedicated core. All other
	// remaining threads are scheduled on the remaining cores using
	// either G-EDF or P-EDF.
	if (algo_ == SchedulingAlgo::GEDF) 
		return schedule_threads_gedf(remain_cores);
	else 
		return schedule_threads_pedf(remain_cores);
}

// Test whether the resulting sequential threads is schedulable with G-EDF.
bool StretchingAlgo::schedule_threads_gedf(uint_t m) {
	//cout << "No. cores for sequential threads: " << m << endl;
	double delta_sum = 0.0, delta_max = 0.0;

	//	uint_t i = 0;
	for (const vector<SequentialThread>& tmp : threads_) {
		for (const SequentialThread& t : tmp) {
			delta_sum += t.work/t.deadline;
			delta_max = max(delta_max, t.work/t.deadline);
		}
		//cout << "i = " << ++i << " -- delta_sum: " << delta_sum << ". delta_max: " << delta_max << endl;
	}

	//	cout << "Final -- delta_sum: " << delta_sum << ". delta_max: " << delta_max << endl;
	if (delta_sum <= (m - (m-1)*delta_max))
		return true;
	return false;
}

// Test whether the resulting sequential threads is schedulable with P-EDF.
// Partitioning is performed with Best Fit Decreasing packing and
// an approximation of demand bound function is used to do uniprocessor
// schedulability test.
bool StretchingAlgo::schedule_threads_pedf(uint_t m) {
	vector<SequentialThread> threads;
	for (const vector<SequentialThread>& v : threads_) {
		threads.insert(threads.end(), v.begin(), v.end());
	}
	
	priority_queue<StretchingPartition, vector<StretchingPartition>, StretchingBestFit> procs;
	for (int i = 0; i < m; ++i) {
		procs.push(StretchingPartition());
	}
	// Process the threads in non-decreasing order of their deadlines.
	sort(threads.begin(), threads.end(), StretchingNonDecreasingDeadline());

	for (int i = 0; i < threads.size(); ++i) {
		if (procs.empty())
			return false;

		bool success = false;
		// Temporarily store partitions that do not fit.
		vector<StretchingPartition> temp;
		while (!procs.empty()) {
			StretchingPartition top_par = procs.top();
			procs.pop();
			if (dbf_approx_partition(top_par, threads[i]) == false) {
				temp.push_back(top_par);
			} else {
				procs.push(top_par);
				for (StretchingPartition& p : temp) {
					procs.push(p);
				}
				success = true;
				break;
			}
		}
		
		if (success == false) {
			for (StretchingPartition& p : temp) {
				procs.push(p);
			}
			return false;
		}
	}
	return true;
}

// Return true if the thread can be inserted to the given processor.
bool StretchingAlgo::dbf_approx_partition(StretchingPartition& par, const SequentialThread& task) {
	// Total demand of all tasks already assigned to this processor.
	double total_dbf = 0;
	for (SequentialThread& t : par.threads) {
		double tmp = Common::dbf_approx(t.work, t.deadline, t.period, task.deadline);
		total_dbf += tmp;
	}

	if (task.deadline >= total_dbf && task.deadline - total_dbf >= task.work) {
		par.threads.push_back(task);
		par.total_load += task.work/task.period;
		return true;
	}
	return false;
}

// Print the MTS representation of the first heavy task.
void StretchingAlgo::print_mts() const {
	const mts_t& mts = mtss_[0];
	for (uint_t i = 0; i < mts.size(); ++i) {
		cout << "Segment " << i << " -- Start: " << mts[i].start
			 << ". Length: " << mts[i].len
			 << ". Subtasks: ";
		for (uint_t node_id : mts[i].subtasks) {
			cout << node_id << " ";
		}
		cout << endl;
	}
}

// Print the resulting threads of the first heavy task.
void StretchingAlgo::print_stretched_threads() const {
	const vector<SequentialThread>& threads = threads_[0];
	for (uint_t i = 0; i < threads.size(); ++i) {
		cout << "Thread " << i << " -- Start: " << threads[i].offset
			 << ". Work: " << threads[i].work << ". Deadline: " << threads[i].deadline << endl;
	}
}
