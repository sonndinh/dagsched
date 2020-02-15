#include "reservation_federated.h"

ReservationFederated::ReservationFederated(const Taskset& taskset, ReservationAlgo reser, Heuristic h, UniScheduler sched, EdfSchedTest test) :
	tset_{taskset}, reser_algo_{reser}, heuristic_{h}, sched_algo_{sched}, sched_test_{test} {
	uint_t no_procs = tset_.get_no_procs();
	procs_.resize(no_procs);
	for (uint_t i = 0; i < no_procs; ++i) {
		procs_[i].id = i; // Proc ID starts from 0.
		procs_[i].load = 0;
	}

	// Initialize priority queue.
	if (heuristic_ == Heuristic::WORST_FIT) {
		for (uint_t i = 0; i < no_procs; ++i) {
			procs_wf_.push(&procs_[i]);
		}
	} else if (heuristic_ == Heuristic::BEST_FIT) {
		for (uint_t i = 0; i < no_procs; ++i) {
			procs_bf_.push(&procs_[i]);
		}
	}
}

// Schedulability test with the given options. 
bool ReservationFederated::is_schedulable() {
	reservation_assign();
	return split_on_fail();
}

// Assign reservation servers to tasks using R-MIN method.
void ReservationFederated::reservation_assign() {
	const vector<DagTask>& heavy_tasks = tset_.get_heavy_tasks();
	const vector<SequentialTask>& light_tasks = tset_.get_light_tasks();
	servers_.resize(heavy_tasks.size() + light_tasks.size());

	uint_t idx = 0;
	if (reser_algo_ == ReservationAlgo::R_MIN) {
		// Compute the servers for heavy tasks.
		for (uint_t i = 0; i < heavy_tasks.size(); ++i, ++idx) {
			uint_t work = heavy_tasks[i].get_work();
			uint_t span = heavy_tasks[i].get_span();
			uint_t deadline = heavy_tasks[i].get_deadline();
			uint_t m_i = ceil(((double)(work - span))/(deadline - span));
			double gamma_i = 1 + (double)(work - span)/(m_i * span);
			
			servers_[idx].deadline = deadline;
			servers_[idx].period = heavy_tasks[i].get_period();
			servers_[idx].work = work;
			servers_[idx].span = span;
			servers_[idx].no_servers = m_i;
			servers_[idx].bound = 2 * m_i;
			servers_[idx].budget = ceil(gamma_i * span);
		}
		
		// Compute the servers for light tasks, one for each light task.
		for (uint_t i = 0; i < light_tasks.size(); ++i, ++idx) {
			servers_[idx].deadline = light_tasks[i].get_deadline();
			servers_[idx].period = light_tasks[i].get_period();
			servers_[idx].work = light_tasks[i].get_work();
			servers_[idx].span = light_tasks[i].get_span();
			servers_[idx].no_servers = 1;
			servers_[idx].bound = 1;
			servers_[idx].budget = light_tasks[i].get_work();
		}
	} else { // ReservationAlgo::R_EQUAL
		// A common stretch ratio for all tasks.
		double gamma = numeric_limits<double>::max();
		for (uint_t i = 0; i < heavy_tasks.size(); ++i) {
			uint_t deadline = heavy_tasks[i].get_deadline();
			uint_t span = heavy_tasks[i].get_span();
			if ((double)deadline/span != 1.0) { // Ignore ratios of 1.0.
				gamma = min(gamma, (double)deadline/span);
			}
		}
		for (uint_t i = 0; i < light_tasks.size(); ++i) {
			uint_t deadline = light_tasks[i].get_deadline();
			uint_t span = light_tasks[i].get_span();
			if ((double)deadline/span != 1.0) { // Ignore ratios of 1.0.
				gamma = min(gamma, (double)deadline/span);
			}
		}

		// If gamma does not change, meaning that all tasks have stretch ratios of 1.0,
		// we assign each task a no. of servers equal to twice its average parallelism, i.e., work/span.
		// This case should happens very rarely.
		if (gamma == numeric_limits<double>::max()) {
			uint_t idx = 0;
			for (uint_t i = 0; i < heavy_tasks.size(); ++i, ++idx) {
				add_requal_servers_special_case(&heavy_tasks[i], idx);
			}
			for (uint_t i = 0; i < light_tasks.size(); ++i, ++idx) {
				add_requal_servers_special_case(&light_tasks[i], idx);
			}
			return;
		}

		// Otherwise, gamma > 1.0.
		uint_t idx = 0;
		// Compute reservation servers for the tasks.
		for (uint_t i = 0; i < heavy_tasks.size(); ++i, ++idx) {
			add_requal_servers(&heavy_tasks[i], idx, gamma);
		}
		for (uint_t i = 0; i < light_tasks.size(); ++i, ++idx) {
			add_requal_servers(&light_tasks[i], idx, gamma);
		}
	}
}

// Adding servers for a DAG task with R-EQUAL algorithm when gamma is 1.0. 
void ReservationFederated::add_requal_servers_special_case(const Task *task, uint_t idx) {
	uint_t work = task->get_work();
	uint_t span = task->get_span();
	uint_t deadline = task->get_deadline();

	servers_[idx].deadline = deadline;
	servers_[idx].period = task->get_period();
	servers_[idx].work = work;
	servers_[idx].span = span;
	uint_t m_i = 2 * ceil(work/span);
	servers_[idx].no_servers = m_i;
	servers_[idx].bound = 2 * m_i;
	servers_[idx].budget = ceil((double)(work + (m_i - 1)*span)/m_i);
}

void ReservationFederated::add_requal_servers(const Task *task, uint_t idx, double gamma) {
	uint_t work = task->get_work();
	uint_t span = task->get_span();
	uint_t deadline = task->get_deadline();
	
	servers_[idx].deadline = deadline;
	servers_[idx].period = task->get_period();
	servers_[idx].work = work;
	servers_[idx].span = span;
	if (work > (uint_t)(gamma * span)) {
		// Classified as heavy task.
		uint_t m_i = ceil((double)(work - span)/(span * (gamma - 1)));
		servers_[idx].no_servers = m_i;
		servers_[idx].bound = 2 * m_i;
		servers_[idx].budget = ceil(gamma * span);
	} else { // Classified as light task.
		servers_[idx].no_servers = 1;
		servers_[idx].bound = 1;
		servers_[idx].budget = work;
	}
}

// Implement the split-on-fail algorithm in the paper (Algorithm 2). 
bool ReservationFederated::split_on_fail() {
	// Sort in non-decreasing order of relative deadline, independent of 
	// whether worst-fit, best-fit, or first-fit heuristic is used.
	sort(servers_.begin(), servers_.end(), DMOrder());

	for (uint_t i = 0; i < servers_.size(); ++i) {
		TaskServers& servers = servers_[i];
		if (servers.no_servers == 1) { // Light task.
			int ret = partition_single(servers.budget, servers.deadline, servers.period);
			if (ret == -1)
				return false;
			
			// Succeeded! Update its processor map.
			servers.server_to_proc.insert({0, ret});
		} else { // Heavy task.
			while (partition_all(servers) == false) {
				if (servers.no_servers > servers.bound)
					return false;
			}
		}
	}
	
	return true;
}

// Partition all reservation servers of a heavy task, starting from the first reservation.
// Return false if it fails to partition any reservation.
// Return true if it partition all reservations successfully.
bool ReservationFederated::partition_all(TaskServers &servers) {
	for (uint_t j = 0; j < servers.no_servers; ++j) {
		int ret = partition_single(servers.budget, servers.deadline, servers.period);
		if (ret == -1) {
			// Revoke all already scheduled servers for this task.
			for (pair<uint_t, uint_t> p : servers.server_to_proc) {
				// Remove the server from the processor.
				procs_[p.second].servers_in_proc.pop_back();
				if (sched_test_ == EdfSchedTest::DBF_APPROX) {
					procs_[p.second].load -= (double)servers.budget/servers.period;
				} else { // EdfSchedTest::DENSITY
					procs_[p.second].load -= (double)servers.budget/servers.deadline;
				}
			}

			// Increase the number of reservations by 1 and update budgets.
			++servers.no_servers;
			servers.budget = ceil((double)(servers.work + (servers.no_servers - 1)*servers.span)/servers.no_servers);
			
			// Clear the map from servers to processors.
			servers.server_to_proc.clear();
			return false;
		}
		
		// Partition succeeded.
		servers.server_to_proc.insert({j, ret});
	}
	
	return true;
}

// Partition a single reservation server with a given heuristic.
// We only support EDF for uniprocessor scheduling for now. 
// Return ID of processor to which the reservation is partitioned to.
// Return -1 if it fails.
int ReservationFederated::partition_single(uint_t budget, uint_t deadline, uint_t period) {
	if (heuristic_ == Heuristic::WORST_FIT) {
		return worst_fit_partition(budget, deadline, period);
	} else if (heuristic_ == Heuristic::BEST_FIT) {
		return best_fit_partition(budget, deadline, period);
	} else { // Heuristic::FIRST_FIT
		return first_fit_partition(budget, deadline, period);
	}
}

// Partition the reservation using worst-fit heuristic.
int ReservationFederated::worst_fit_partition(uint_t budget, uint_t deadline, uint_t period) {
	ReserFedProcessor *top = procs_wf_.top();
	procs_wf_.pop();
	if (can_partition(top, budget, deadline, period) == false) {
		procs_wf_.push(top);
		return -1;
	}

	// Partition succeeded!
	procs_wf_.push(top);
	return top->id;
}

// Partition the reservation using best-fit heuristic.
int ReservationFederated::best_fit_partition(uint_t budget, uint_t deadline, uint_t period) {
	// Temporarily store processors that do not have enough room.
	vector<ReserFedProcessor*> tmp;
	while (!procs_bf_.empty()) {
		ReserFedProcessor *top = procs_bf_.top();
		procs_bf_.pop();
		
		if (can_partition(top, budget, deadline, period) == true) {
			procs_bf_.push(top);
			// Insert back the temporarily stored processors.
			for (ReserFedProcessor *p : tmp) {
				procs_bf_.push(p);
			}
			return top->id;
		} else {
			tmp.push_back(top);
		}
	}

	// Fail to partition this reservation.
	for (ReserFedProcessor *p : tmp) {
		procs_bf_.push(p);
	}
	return -1;
}

// Partition the reservation using first-fit heuristic.
int ReservationFederated::first_fit_partition(uint_t budget, uint_t deadline, uint_t period) {
	for (uint_t i = 0; i < procs_.size(); ++i) {
		if (can_partition(&procs_[i], budget, deadline, period) == true) {
			return procs_[i].id;
		}
	}
	return -1;
}

// Run schedulability test on a processor with a new reservation.
// Return true if the reservation can be added to this processor. 
// Return false otherwise.
bool ReservationFederated::can_partition(ReserFedProcessor *proc, uint_t budget, uint_t deadline, uint_t period) {
	if (sched_test_ == EdfSchedTest::DENSITY) {
		// Use DENSITY test.
		if (proc->load + (double)budget/deadline <= 1.0) {
			proc->servers_in_proc.push_back({budget, deadline, period});
			proc->load += (double)budget/deadline;
			return true;
		}
		return false;
	}

	// Otherwise, use EdfSchedTest::DBF_APPROX test.
	// Total demand of all servers already assigned to this processor.
	uint_t total_dbf = 0;
	for (Server s : proc->servers_in_proc) {
		uint_t tmp = Common::dbf_approx(s.budget, s.deadline, s.period, deadline);
		total_dbf += tmp;
	}

	if (deadline >= total_dbf && deadline - total_dbf >= budget) {
		// Insert this server to the processor.
		proc->servers_in_proc.push_back({budget, deadline, period});
		proc->load += (double)budget/period;
		return true;
	}
	
	return false;
}
