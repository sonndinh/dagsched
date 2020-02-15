#include <iostream>
#include <queue>
#include "semi_federated.h"

SemiFederated::SemiFederated(const Taskset& taskset) : tset_{taskset} {
	avail_procs_ = tset_.get_no_procs();
	heavy_procs_ = 0;
}


// Implement the SF[x+2] algorithm in the semi-federated paper.
bool SemiFederated::is_schedulable() {
    if (schedule_heavy_tasks() == false)
		return false;

	if (schedule_seq_tasks() == false)
		return false;
	return true;
}

uint_t SemiFederated::procs_to_heavy_tasks() const {
	return heavy_procs_;
}

// Return FALSE if there are not enough cores to allocate for heavy tasks.
// Return TRUE otherwise.
bool SemiFederated::schedule_heavy_tasks() {
	const vector<DagTask>& heavy_tasks = tset_.get_heavy_tasks();
	
	for (int i = 0; i < heavy_tasks.size(); ++i) {
		uint_t work = heavy_tasks[i].get_work();
		uint_t span = heavy_tasks[i].get_span();
		uint_t deadline = heavy_tasks[i].get_deadline();
		if (deadline == span) {
			// Need infinite number of cores.
			return false;
		} else {
			// Heavy tasks have density > 1.0. Thus work > deadline.
			// Since deadline >= span, work > span.
			double required_cap = (double)(work - span)/(deadline - span);
			uint_t dedicated_cores = (work - span)/(deadline - span);
			//			cout << "[Semi-Federated] Heavy task " << i << " needs " << dedicated_cores << endl;
			if (avail_procs_ < dedicated_cores) {
				return false;
			}
			avail_procs_ -= dedicated_cores;
			heavy_procs_ += dedicated_cores;

			// Create container task for this heavy task.
			double load_bound = required_cap - (double)dedicated_cores;
			double load_star = max(load_bound/2, load_bound/required_cap);
			container_tasks_.push_back(SequentialTask(load_bound, load_star));
		}
	}
	
	return true;
}

// Return FALSE if the remaining cores cannot schedule light tasks and
// container tasks. Return TRUE otherwise.
// This implements Algorithm 3 in the paper.
bool SemiFederated::schedule_seq_tasks() {
	vector<SequentialTask> seq_tasks(tset_.get_light_tasks());
	for (int i = 0; i < container_tasks_.size(); ++i) {
		seq_tasks.push_back(container_tasks_[i]);
	}
	// Sort all sequential tasks, including light tasks and container tasks
	// in non-increasing order of their delta* loads.
	sort(seq_tasks.begin(), seq_tasks.end(), NonIncreasingLoadStar());

	// List of full partition which cannot receive any further tasks.
	vector<Partition> full_procs;
	// Priority queue of available partitions. The one with smallest
	// total delta* load is on top, i.e., for Worst-Fit packing.
	priority_queue<Partition, vector<Partition>, SmallerLoadStar> avail_procs;
	for (int i = 0; i < avail_procs_; ++i) {
		avail_procs.push(Partition());
	}

	// Pack the sequential tasks in worst-fit order.
	for (int i = 0; i < seq_tasks.size(); ++i) {
		if (avail_procs.empty())
			return false;
		double top_load = avail_procs.top().total_load_star;
		if (top_load + seq_tasks[i].get_load_star() > 1.0)
			return false;
		
		// Put the task to this processor.
		Partition top = avail_procs.top();
		avail_procs.pop();
		top.tasks.push_back(seq_tasks[i]);
		top.total_load += seq_tasks[i].get_load();
		top.total_load_star += seq_tasks[i].get_load_star();

		if (top.total_load > 1.0) {
			// This processor is full, move it to the full list.
			full_procs.push_back(top);
		} else {
			// Otherwise, push back to the priority queue.
			avail_procs.push(top);
		}
	}

	seq_tasks.clear();
	// Scrape out the exceeding loads from the full processors.
	for (int i = 0; i < full_procs.size(); ++i) {
		vector<SequentialTask> temp = scrape(full_procs[i]);
		seq_tasks.insert(seq_tasks.end(), temp.begin(), temp.end());
	}

	// Partition the sequential tasks obtained after scraping into
	// the processors that are not full in worst-fit order.
	sort(seq_tasks.begin(), seq_tasks.end(), NonIncreasingLoad());
	// Sort the available processors in minimal total load, i.e.,
	// processor with smallest total load is on top.
	priority_queue<Partition, vector<Partition>, SmallerLoad> remain_procs;
	while (!avail_procs.empty()) {
		remain_procs.push(avail_procs.top());
		avail_procs.pop();
	}

	// Worst-Fit Packing.
	for (int i = 0; i < seq_tasks.size(); ++i) {
		if (remain_procs.empty())
			return false;
		Partition top = remain_procs.top();
		// Fail if the processor with minimal load cannot handle this task.
		// Note that the remaining load for the scraped task is store in
		// its @load_ attribute.
		if (top.total_load + seq_tasks[i].get_load() > 1.0)
			return false;

		remain_procs.pop();
		top.tasks.push_back(seq_tasks[i]);
		top.total_load += seq_tasks[i].get_load();
		remain_procs.push(top);
	}
	
	return true;
}

// Scrape exceeding load from a full partition.
// This implements the Algorithm 4 in the paper.
vector<SequentialTask> SemiFederated::scrape(Partition& par) {
	vector<SequentialTask> scraped_loads;
	double exceeding_load = par.total_load - 1.0;

	// Go through the container tasks and see if we can scrape a part of them out.
	for (int i = 0; i < par.tasks.size(); ++i) {
		SequentialTask &task = par.tasks[i];
		if (task.get_type() == SequentialTask::SeqTaskType::LIGHT)
			continue;

		if (task.get_load() - task.get_load_star() > exceeding_load) {
			SequentialTask remain = task;
			remain.set_load(exceeding_load);
			task.set_load(task.get_load() - exceeding_load);
			scraped_loads.push_back(remain);
			par.total_load = 1.0;
			return scraped_loads;
		} else {
			SequentialTask remain = task;
			double diff = task.get_load() - task.get_load_star();
			remain.set_load(diff);
			task.set_load(task.get_load_star());
			scraped_loads.push_back(remain);
			par.total_load -= diff;
			exceeding_load -= diff;
		}
	}

	if (par.total_load > 1.0) {
		cout << "WARNING: A processor still has load > 1.0 after scraping !!" << endl;
	}
	return scraped_loads;
}
