#include <queue>
#include "federated_common.h"

FederatedScheduler::FederatedScheduler(const Taskset& taskset, Heuristic h, EdfSchedTest test)
	: tset_{taskset}, heuristic_{h}, sched_test_{test} {
	avail_procs_ = taskset.get_no_procs();
	heavy_procs_ = 0;
}

bool FederatedScheduler::is_schedulable() {
	if (schedule_heavy_tasks() == false)
		return false;
	if (schedule_light_tasks() == false)
		return false;
	return true;
}

bool FederatedScheduler::schedule_light_tasks() {
	if (heuristic_ == Heuristic::WORST_FIT) {
		return worst_fit_partition();
	} else if (heuristic_ == Heuristic::BEST_FIT) {
		return best_fit_partition();
	} else { // Heuristic::FIRST_FIT
		return first_fit_partition();
	}	
}

// Partition in Worst-Fit order.
bool FederatedScheduler::worst_fit_partition() {
	vector<SequentialTask> light_tasks = tset_.get_light_tasks();
	priority_queue<ProposedPartition, vector<ProposedPartition>, WorstFitCompare> procs;
	for (int i = 0; i < avail_procs_; ++i) {
		procs.push(ProposedPartition());
	}

	// Sort the light tasks in non-decreasing order of deadline if use
	// DBF* function for EDF schedulability test.
	// Sort them using their densities if use the total density test.
	if (sched_test_ == EdfSchedTest::DBF_APPROX) {
		sort(light_tasks.begin(), light_tasks.end(), NonDecreasingDeadline());
	} else { // EdfSchedTest::DENSITY
		sort(light_tasks.begin(), light_tasks.end(), SortWorstFit());
	}

	// Partition each light task.
	for (int i = 0; i < light_tasks.size(); ++i) {
		if (procs.empty())
			return false;
		
		ProposedPartition top_par = procs.top();
		if (sched_test_ == EdfSchedTest::DBF_APPROX) {
			if (dbf_approx_partition(top_par, light_tasks[i]) == false)
				return false;
		} else { // EdfSchedTest::DENSITY
			// Check if total density > 1.0.
			if (top_par.total_load + light_tasks[i].get_load() > 1.0)
				return false;
			top_par.total_load += light_tasks[i].get_load();
			top_par.tasks.push_back(light_tasks[i]);
		}
		
		procs.pop();
		procs.push(top_par);
	}
	
	return true;
}

// Partition in Best-Fit order.
bool FederatedScheduler::best_fit_partition() {
	vector<SequentialTask> light_tasks = tset_.get_light_tasks();
	priority_queue<ProposedPartition, vector<ProposedPartition>, BestFitCompare> procs;
	for (int i = 0; i < avail_procs_; ++i) {
		procs.push(ProposedPartition());
	}
	
	if (sched_test_ == EdfSchedTest::DBF_APPROX) {
		sort(light_tasks.begin(), light_tasks.end(), NonDecreasingDeadline());
	} else { // EdfSchedTest::DENSITY
		sort(light_tasks.begin(), light_tasks.end(), SortBestFit());
	}

	for (int i = 0; i < light_tasks.size(); ++i) {
		if (procs.empty())
			return false;

		bool success = false;
		// Temporarily store partitions that do not fit.
		vector<ProposedPartition> temp;
		while (!procs.empty()) {
			ProposedPartition top_par = procs.top();
			procs.pop();
			if (sched_test_ == EdfSchedTest::DBF_APPROX) {
				if (dbf_approx_partition(top_par, light_tasks[i]) == false) {
					temp.push_back(top_par);
				} else {
					procs.push(top_par);
					for (ProposedPartition& p : temp) {
						procs.push(p);
					}
					success = true;
					break;
				}
			} else { // EdfSchedTest::DENSITY
				if (top_par.total_load + light_tasks[i].get_load() > 1.0) {
					temp.push_back(top_par);
				} else {
					top_par.total_load += light_tasks[i].get_load();
					top_par.tasks.push_back(light_tasks[i]);
					procs.push(top_par);
					for (ProposedPartition& p : temp) {
						procs.push(p);
					}
					success = true;
					break;
				}
			}
		}
		
		if (success == false) {
			for (ProposedPartition& p : temp) {
				procs.push(p);
			}
			return false;
		}
	}
	return true;
}

// Partition in First-Fit order.
bool FederatedScheduler::first_fit_partition() {
	vector<SequentialTask> light_tasks = tset_.get_light_tasks();
	sort(light_tasks.begin(), light_tasks.end(), NonDecreasingDeadline());
	vector<ProposedPartition> procs(avail_procs_);

	for (int i = 0; i < light_tasks.size(); ++i) {
		bool success = false;
		for (int j = 0; j < procs.size(); ++j) {
			if (sched_test_ == EdfSchedTest::DBF_APPROX) {
				if (dbf_approx_partition(procs[j], light_tasks[i]) == true) {
					success = true;
					break;
				}
			} else { // EdfSchedTest::DENSITY
				if (procs[j].total_load + light_tasks[i].get_load() <= 1.0) {
					procs[j].total_load += light_tasks[i].get_load();
					procs[j].tasks.push_back(light_tasks[i]);
					success = true;
					break;
				}
			}
		}
		if (success == false) {
			return false;
		}
	}

	return true;
}

// Attempt to partition a light task in to a processor.
// Return true if succeeds. Return false if fails.
bool FederatedScheduler::dbf_approx_partition(ProposedPartition& par, const SequentialTask& task) {
	// Total demand of all tasks already assigned to this processor.
	double total_dbf = 0;
	for (SequentialTask& t : par.tasks) {
		double tmp = Common::dbf_approx(t.get_work(), t.get_deadline(), t.get_period(), task.get_deadline());
		total_dbf += tmp;
	}

	if ((double)task.get_deadline() >= total_dbf && (double)task.get_deadline() - total_dbf >= (double)task.get_work()) {
		par.tasks.push_back(task);
		par.total_load += task.get_utilization();
		return true;
	}
	return false;
}
