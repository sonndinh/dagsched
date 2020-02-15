// Implementation for the algorithm in paper "Reservation-Based Federated
// Scheduling for Parallel Real-Time Tasks" by Ueter et al., published at RTSS 2018.
#ifndef __RESERVATION_FED
#define __RESERVATION_FED

#include <vector>
#include <unordered_map>
#include <queue>
#include "taskset.h"
#include "common.h"

// Scheduling algorithm for partitions.
// Only support EDF for now.
enum class UniScheduler {EDF};

// Two algorithms to assign reservations.
enum class ReservationAlgo {R_MIN, R_EQUAL};

// Store information for reservation servers of each DAG task.
// Light tasks receive 1 reservation server each.
// Heavy tasks receive > 1 reservation servers each.
struct TaskServers {
	// Deadline and period of the corresponding DAG task.
	uint_t deadline;
	uint_t period;
	// Work and span of the corresponding DAG task.
	// These are used to recompute the budget.
	uint_t work;
	uint_t span;
	
	// Number of reservation servers for this task.
	uint_t no_servers;
	// Boundary number for the Split-On-Fail algorithm.
	uint_t bound;
	
	// Budget of each reservation server.
	uint_t budget;

	// Map from reservation to the processor it is partitioned to.
	unordered_map<uint_t, uint_t> server_to_proc;
};

// For sorting the reservation servers based on deadline-monotonic order.
struct DMOrder {
	bool operator()(const TaskServers& a, const TaskServers& b) {
		return a.deadline < b.deadline;
	}
};

// A simple representation for reservation server for mapping to processor.
struct Server {
	uint_t budget;
	uint_t deadline;
	uint_t period;
};

// Representing a processor to which reservations are partitioned.
struct ReserFedProcessor {
	uint_t id;
	
	// For the schedulability test that uses DBF-approx function,
	// this is the total utilizations of all servers assigned to this processor.
	// For schedulability test using total density, this is the total
	// density of all servers already assigned to this processor.
	double load;
	vector<Server> servers_in_proc;
};

// Processors with smaller total loads get popped first.
struct WorstFitOrder {
	bool operator()(const ReserFedProcessor *a, const ReserFedProcessor *b) {
		return a->load > b->load;
	}
};

// Processors with larger total loads get popped first.
struct BestFitOrder {
	bool operator()(const ReserFedProcessor *a, const ReserFedProcessor *b) {
		return a->load < b->load;
	}
};

class ReservationFederated {
public:
	ReservationFederated(const Taskset& taskset, ReservationAlgo reser, Heuristic h,
						 UniScheduler sched = UniScheduler::EDF, EdfSchedTest test = EdfSchedTest::DBF_APPROX);
	bool is_schedulable();

private:
	// Convert DAG tasks to reservation servers.
	void reservation_assign();
	void add_requal_servers(const Task* task, uint_t idx, double gamma);
	void add_requal_servers_special_case(const Task* task, uint_t idx);

	// Partition a single reservation server to a processor.
	int partition_single(uint_t budget, uint_t deadline, uint_t period);
	int worst_fit_partition(uint_t, uint_t, uint_t);
	int best_fit_partition(uint_t, uint_t, uint_t);
	int first_fit_partition(uint_t, uint_t, uint_t);
	bool can_partition(ReserFedProcessor*, uint_t, uint_t, uint_t);

	// Partition all reservation servers belonging to a DAG task.
	bool partition_all(TaskServers &servers);

	// The main function for the Split-On-Fail function.
	bool split_on_fail();
	
private:
	Taskset tset_;
	ReservationAlgo reser_algo_;
	Heuristic heuristic_;
	UniScheduler sched_algo_;
	EdfSchedTest sched_test_;
	
	// Reservation servers obtained after calling R-MIN or R-EQUAL.
	vector<TaskServers> servers_;

	// Processors of the system.
	vector<ReserFedProcessor> procs_;

	// Queue of processors prioritized for worst-fit or best-fit.
	priority_queue<ReserFedProcessor*, vector<ReserFedProcessor*>, WorstFitOrder> procs_wf_;
	priority_queue<ReserFedProcessor*, vector<ReserFedProcessor*>, BestFitOrder> procs_bf_;
};

#endif
