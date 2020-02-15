#include <iostream>
#include <fstream>
#include <cassert>
#include "taskset.h"

//void print_dag_task(const DagTask&);

// @m: Total number of processors in the system.
// @normalized_util: Total normalized utilization.
// @heavy_ratio: Percentage of total utilization occupied by heavy tasks.
// @num_heavy: Number of heavy tasks.
// @heavy_util_file: Path to a file containing individual utilizations for heavy tasks.
// @size_min, @size_max: Generate numbers of nodes for heavy tasks in [size_min, size_max].
// @wcet_min, @wcet_max: Generate wcets of nodes for heavy tasks in [wcet_min, wcet_max].
// @p: Probability for an edge between any 2 nodes.
// @num_light: Number of light tasks.
// @light_util_file: Path to a file containing utilizations for light tasks.
// NOTE: We consider a task heavy if it has density > 1.0, light if it has density <= 1.0.
Taskset::Taskset(uint_t m, double normalized_util, double heavy_ratio, uint_t num_heavy,
				 string heavy_util_file, uint_t size_min, uint_t size_max, uint_t wcet_min,
				 uint_t wcet_max, double p, uint_t num_light, string light_util_file) : num_procs_(m) {
	total_util_ = m * normalized_util;
	heavy_util_ = heavy_ratio * total_util_;

	// Read utilizations for heavy tasks from file.
	vector<double> heavy_utils;
	fstream heavy_f(heavy_util_file, std::fstream::in);
	string line;
	double temp;
	while (getline(heavy_f, line)) {
		temp = stod(line);
		heavy_utils.push_back(temp);
	}
	assert(heavy_utils.size() == num_heavy);
	heavy_f.close();

	// Generate heavy tasks.
	heavy_tasks_.clear();
	for (int i = 0; i < num_heavy; ++i) {
		DagTask task(size_min, size_max, wcet_min, wcet_max, p, heavy_utils[i]);
		//		print_dag_task(task);
		heavy_tasks_.push_back(task);
	}
	
	// Read utilizations for light tasks from file and generate.
	vector<double> light_utils;
	fstream light_f(light_util_file, std::fstream::in);
	while (getline(light_f, line)) {
		temp = stod(line);
		light_utils.push_back(temp);
	}
	assert(light_utils.size() == num_light);
	light_f.close();

	// Generate light tasks.
	light_tasks_.clear();
	for (int i = 0; i < num_light; ++i) {
		light_tasks_.push_back(SequentialTask(size_min, size_max, wcet_min, wcet_max, p, light_utils[i]));
	}
}


// Generate tasks with given densities instead of utilizations.
// @dummy: Just an additional argument to distinguish with the other constructor.
Taskset::Taskset(uint_t m, double normalized_den, double heavy_ratio, uint_t num_heavy,
				 string heavy_file, uint_t size_min, uint_t size_max, uint_t wcet_min,
				 uint_t wcet_max, double p, uint_t num_light, string light_file, bool dummy) : num_procs_(m) {
	// We don't care about utilization for this case.
	total_util_ = 0;
	heavy_util_ = 0;

	// Read densities for heavy tasks from file.
	vector<double> heavy_densities;
	fstream heavy_f(heavy_file, std::fstream::in);
	string line;
	double temp;
	while (getline(heavy_f, line)) {
		temp = stod(line);
		heavy_densities.push_back(temp);
	}
	//	cout << "heavy_densities.size(): " << heavy_densities.size() << ". num_heavy: " << num_heavy << endl;
	assert(heavy_densities.size() == num_heavy);
	heavy_f.close();

	// Generate heavy tasks.
	heavy_tasks_.clear();
	for (int i = 0; i < num_heavy; ++i) {
		DagTask task(size_min, size_max, wcet_min, wcet_max, p, heavy_densities[i], true);
		//		print_dag_task(task);
		heavy_tasks_.push_back(task);
	}
	
	// Read utilizations for light tasks from file and generate.
	vector<double> light_densities;
	fstream light_f(light_file, std::fstream::in);
	while (getline(light_f, line)) {
		temp = stod(line);
		light_densities.push_back(temp);
	}
	assert(light_densities.size() == num_light);
	light_f.close();

	// Generate light tasks.
	light_tasks_.clear();
	for (int i = 0; i < num_light; ++i) {
		light_tasks_.push_back(SequentialTask(size_min, size_max, wcet_min, wcet_max, p, light_densities[i], true));
	}
}

// Change the range for D_i/T_i for all tasks in this task set.
// Consequently, each task will generate a new period.
void Taskset::change_deadline_to_period_ratio(double min, double max) {
	for (int i = 0; i < heavy_tasks_.size(); ++i) {
		heavy_tasks_[i].change_deadline_to_period_ratio(min, max);
	}

	for (int i = 0; i < light_tasks_.size(); ++i) {
		light_tasks_[i].change_deadline_to_period_ratio(min, max);
	}
}
