#ifndef __TASK
#define __TASK

#include "common.h"

// Base class for classes DagTask and SequentialTask.
class Task {
public:
	Task(uint_t work = 0, uint_t span = 0, uint_t period = 0, uint_t deadline = 0) :
		work_{work}, span_{span}, deadline_{deadline}, period_{period} {}
	
	uint_t get_work() const {
		return work_;
	}

	uint_t get_span() const {
		return span_;
	}

	uint_t get_deadline() const {
		return deadline_;
	}

	uint_t get_period() const {
		return period_;
	}

	double get_utilization() const {
		return (double)work_/period_;
	}

	virtual uint_t get_size() const = 0;
	virtual const matrix_t& get_adj_matrix() const = 0;
	virtual const std::vector<uint_t>& get_node_lengths() const = 0;
	virtual const std::vector<uint_t>& get_crit_path() const = 0;
	
protected:
	uint_t work_;
	uint_t span_;
	uint_t deadline_;
	uint_t period_;
};

#endif
