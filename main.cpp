// A bunch of main functions and test cases for evaluation.
#include <vector>
#include <iostream>
#include <string>
#include <stack>
#include <fstream>
#include <sstream>
#include <cassert>
#include "proposed_scheduler.h"
#include "naive_greedy_scheduler.h"
#include "baruah_greedy_scheduler.h"
#include "taskset.h"
#include "semi_federated.h"
#include "reservation_federated.h"
#include "baruah_federated.h"
#include "proposed_federated.h"
#include "stretching_algorithm.h"
#include "decomp_algorithm.h"
using namespace std;


void print_proposed_schedule(const ProposedScheduler&);
void print_dag_task(const DagTask&);
void print_naive_schedule(const NaiveGreedyScheduler&, const DagTask&);
void print_baruah_schedule(const BaruahGreedyScheduler&, const DagTask&);
void read_dag_from_file(string file_path, vector<uint_t>& wcets, matrix_t &adj_matrix, uint_t &period, uint_t &deadline);
int test_extra_preemptions(int, char**);
int expr_vary_period_ratio_fullset(int argc, char *argv[]);
int expr_vary_util_fullset(int argc, char *argv[]);

int main(int argc, char* argv[]) {
	//return expr_vary_util_fullset(argc, argv);
	return expr_vary_period_ratio_fullset(argc, argv);
}

// MAIN FUNCTION FOR FULL TASK SET EXPERIMENTS THAT VARY
// THE RATIO OF D_i/T_i.
// Generate task sets of fixed total density using RandFixedSum.
// Control the total heavy task density (and thus the total light task density),
// and the number of heavy tasks (and the number of light tasks).
// For each task, D_i is derived from the generated density and its DAG's C_i.
// Make sure the computed D_i >= L_i.
int expr_vary_period_ratio_fullset(int argc, char *argv[]) {
	if (argc != 14) {
		cout << "Usage: Program <NumProcs> <TotalNormalizedDensity> <HeavyRatio> <NumHeavy> <MinNumNodes> \
<MaxNumNodes> <MinWcet> <MaxWcet> <EdgeProb> <NumLight> <InputPath> <NumTasksets> <OutputPath> !!" << endl;
		return -1;
	}

	// Different range for the ratio D_i/T_i.
	vector<vector<double>> beta_list{/*{0.01, 0.05}, {0.05, 0.1}, {0.1, 0.15}, {0.15, 0.2}, {0.2, 0.25},
									{0.25, 0.3}, {0.3, 0.35}, {0.35, 0.4}, {0.4, 0.45}, {0.45, 0.5},
									{0.5, 0.55},*/ {0.55, 0.6}, {0.6, 0.65}, {0.65, 0.7}, {0.7, 0.75},
									{0.75, 0.8}, {0.8, 0.85}, {0.85, 0.9}, {0.9, 0.95}, {0.95, 1.0}};
	
	uint_t m = stoi(string(argv[1]));
	double normalized_den = stod(string(argv[2]));
	double heavy_ratio = stod(string(argv[3]));
	uint_t num_heavy = stoi(string(argv[4]));
	uint_t size_min = stoi(string(argv[5]));
	uint_t size_max = stoi(string(argv[6]));
	uint_t wcet_min = stoi(string(argv[7]));
	uint_t wcet_max = stoi(string(argv[8]));
	double p = stod(string(argv[9]));
	uint_t num_light = stoi(string(argv[10]));
	uint_t num_tsets = stoi(string(argv[12]));

	// Generate all task sets.
	// The periods for tasks will be updated later depending on the ratio D_i/T_i.
	vector<Taskset> tasksets;
	for (int i = 1; i <= num_tsets; ++i) {
		string tset_path = string(argv[11]) + "/taskset" + to_string(i);
		string heavy_file = tset_path + "/heavy_densities.txt";
		string light_file = tset_path + "/light_densities.txt";
		tasksets.push_back(Taskset(m, normalized_den, heavy_ratio, num_heavy, heavy_file, size_min,
								   size_max, wcet_min, wcet_max, p, num_light, light_file, true));
	}
	
	for (vector<double> beta : beta_list) {
		//uint_t count_semi = 0;
		//uint_t count_baruah_bf_dbf = 0;
		//uint_t count_proposed_wf_den = 0;
		//uint_t count_proposed_bf_dbf = 0;
		//uint_t count_reser_bf_dbf = 0;
		uint_t count_stretching_gedf = 0;
		uint_t count_stretching_pedf = 0;
		uint_t count_decomp_algo = 0;
		cout << "=======> BETA (D_i/T_i) range: [" << beta[0] << ", " << beta[1] << ")" << endl;
		for (int i = 1; i <= num_tsets; ++i) {
			if (i % 50 == 0) {
				cout << "Computing for task set " << i << " ..." << endl;
			}

			// Update the period for the tasks in this task set using the range of beta.
			Taskset &taskset = tasksets[i-1];
			taskset.change_deadline_to_period_ratio(beta[0], beta[1]);

			/*
			SemiFederated semi(taskset);
			BaruahFederated baruah_bf_dbf(taskset, Heuristic::BEST_FIT, EdfSchedTest::DBF_APPROX);
			ProposedFederated proposed_wf_den(taskset, Heuristic::WORST_FIT, EdfSchedTest::DENSITY);
			ProposedFederated proposed_bf_dbf(taskset, Heuristic::BEST_FIT, EdfSchedTest::DBF_APPROX);
			ReservationFederated reser_bf_dbf(taskset, ReservationAlgo::R_MIN, Heuristic::BEST_FIT, UniScheduler::EDF, EdfSchedTest::DBF_APPROX);
			
			if (reser_bf_dbf.is_schedulable()) {
				++count_reser_bf_dbf;
			}
			if (semi.is_schedulable()) {
				++count_semi;
			}
			if (baruah_bf_dbf.is_schedulable()) {
				++count_baruah_bf_dbf;
			}
			if (proposed_wf_den.is_schedulable()) {
				++count_proposed_wf_den;
			}
			if (proposed_bf_dbf.is_schedulable()) {
				++count_proposed_bf_dbf;
			}
			*/
			StretchingAlgo stretch_gedf(taskset, SchedulingAlgo::GEDF);
			StretchingAlgo stretch_pedf(taskset, SchedulingAlgo::PEDF);
			DecompAlgo decomp(taskset);

			if (stretch_gedf.is_schedulable()) {
				++count_stretching_gedf;
			}
			if (stretch_pedf.is_schedulable()) {
				++count_stretching_pedf;
			}
			if (decomp.is_schedulable()) {
				++count_decomp_algo;
			}
		}
		/*
		cout << "Reservation-Based [BF-RMIN-DBF] ==> Acceptance Ratio: " << count_reser_bf_dbf << "/" << num_tsets << endl;
		cout << "Semi-Federated ==> Acceptance Ratio: " << count_semi << "/" << num_tsets << endl;
		cout << "Baruah Algo [BF-DBF] ==> Acceptance Ratio: " << count_baruah_bf_dbf << "/" << num_tsets << endl;
		cout << "Proposed Algo [WF-DEN]  ==> Acceptance Ratio: " << count_proposed_wf_den << "/" << num_tsets << endl;
		cout << "Proposed Algo [BF-DBF]  ==> Acceptance Ratio: " << count_proposed_bf_dbf << "/" << num_tsets << endl;
		
		// Write to a result file.
		fstream out_fs(string(argv[11]) + "/result_beta_" + to_string(beta[0]) + "_" + to_string(beta[1]) + ".txt", fstream::out);
		out_fs << count_reser_bf_dbf << "/" << num_tsets << endl;
		out_fs << count_semi << "/" << num_tsets << endl;
		out_fs << count_baruah_bf_dbf << "/" << num_tsets << endl;
		out_fs << count_proposed_wf_den << "/" << num_tsets << endl;
		out_fs << count_proposed_bf_dbf << "/" << num_tsets << endl;
		out_fs.close();
		*/
		cout << "Stretching Algo [G-EDF]  ==> Acceptance Ratio: " << count_stretching_gedf << "/" << num_tsets << endl;
		cout << "Stretching Algo [P-EDF]  ==> Acceptance Ratio: " << count_stretching_pedf << "/" << num_tsets << endl;
		cout << "Decomposition Algo  ==> Acceptance Ratio: " << count_decomp_algo << "/" << num_tsets << endl;

		fstream out_fs(string(argv[11]) + "/result_stretch_beta_" + to_string(beta[0]) + "_" + to_string(beta[1]) + ".txt", fstream::out);
		out_fs << count_stretching_gedf << "/" << num_tsets << endl;
		out_fs << count_stretching_pedf << "/" << num_tsets << endl;
		out_fs << count_decomp_algo << "/" << num_tsets << endl;
		out_fs.close();
	}
	
	return 0;
}


// MAIN FUNCTION FOR FULL TASK SET EXPERIMENTS.
// Read utilizations for heavy tasks and light tasks from files.
// Generate task sets and run the proposed federated scheduling, semi-federated,
// reservation-based federated scheduling, and baruah's federated scheduling.
int expr_vary_util_fullset(int argc, char *argv[]) {
	if (argc != 14) {
		cout << "Usage: Program <NumProcs> <TotalNormalizedUtil> <HeavyRatio> <NumHeavy> <MinNumNodes> \
<MaxNumNodes> <MinWcet> <MaxWcet> <EdgeProb> <NumLight> <InputPath> <NumTasksets> <OutputPath> !!" << endl;
		return -1;
	}

	uint_t m = stoi(string(argv[1]));
	double normalized_util = stod(string(argv[2]));
	double heavy_ratio = stod(string(argv[3]));
	uint_t num_heavy = stoi(string(argv[4]));	
	uint_t size_min = stoi(string(argv[5]));
	uint_t size_max = stoi(string(argv[6]));
	uint_t wcet_min = stoi(string(argv[7]));
	uint_t wcet_max = stoi(string(argv[8]));
	double p = stod(string(argv[9]));
	uint_t num_light = stoi(string(argv[10]));
	uint_t num_tsets = stoi(string(argv[12]));

	uint_t count_semi = 0;
	uint_t count_baruah_bf_dbf = 0;
	uint_t count_proposed_wf_den = 0;
	//uint_t count_proposed_wf_dbf = 0;
	//uint_t count_proposed_bf_den = 0;
	uint_t count_proposed_bf_dbf = 0;
	uint_t count_reser_bf_dbf = 0;
	//	uint_t count_reser_bf_den = 0;
	//	uint_t count_reser_requal_bf = 0;
	uint_t count_stretching_gedf = 0;
	uint_t count_stretching_pedf = 0;
	uint_t count_decomp_algo = 0;
	for (int i = 1; i <= num_tsets; ++i) {
		if (i % 50 == 0) {
			cout << "Computing for task set " << i << " ..." << endl;
		}
		string tset_path = string(argv[11]) + "/taskset" + to_string(i);
		string heavy_utils_file = tset_path + "/heavy_utils.txt";
		string light_utils_file = tset_path + "/light_utils.txt";
		Taskset taskset(m, normalized_util, heavy_ratio, num_heavy, heavy_utils_file, size_min,
						size_max, wcet_min, wcet_max, p, num_light, light_utils_file);
		
		//SemiFederated semi(taskset);
		//BaruahFederated baruah_bf_dbf(taskset, Heuristic::BEST_FIT, EdfSchedTest::DBF_APPROX);
		//ProposedFederated proposed_wf_den(taskset, Heuristic::WORST_FIT, EdfSchedTest::DENSITY);
		//		ProposedFederated proposed_wf_dbf(taskset, Heuristic::WORST_FIT, EdfSchedTest::DBF_APPROX);
		//		ProposedFederated proposed_bf_den(taskset, Heuristic::BEST_FIT, EdfSchedTest::DENSITY);
		//ProposedFederated proposed_bf_dbf(taskset, Heuristic::BEST_FIT, EdfSchedTest::DBF_APPROX);
		//ReservationFederated reser_bf_dbf(taskset, ReservationAlgo::R_MIN, Heuristic::BEST_FIT, UniScheduler::EDF, EdfSchedTest::DBF_APPROX);
		//		ReservationFederated reser_bf_den(taskset, ReservationAlgo::R_MIN, Heuristic::BEST_FIT, UniScheduler::EDF, EdfSchedTest::DENSITY);		
		//		ReservationFederated reser_requal_bf(taskset, ReservationAlgo::R_EQUAL, Heuristic::BEST_FIT);
		StretchingAlgo stretch_gedf(taskset, SchedulingAlgo::GEDF);
		StretchingAlgo stretch_pedf(taskset, SchedulingAlgo::PEDF);
		DecompAlgo decomp(taskset);

		/*
		if (reser_bf_dbf.is_schedulable()) {
			++count_reser_bf_dbf;
		}
		//if (reser_bf_den.is_schedulable()) {
		//	++count_reser_bf_den;
		//}
		//if (reser_requal_bf.is_schedulable()) {
		//	++count_reser_requal_bf;
		//}
		if (semi.is_schedulable()) {
			++count_semi;
		}
		if (baruah_bf_dbf.is_schedulable()) {
			++count_baruah_bf_dbf;
		}
		if (proposed_wf_den.is_schedulable()) {
			++count_proposed_wf_den;
		}
		//if (proposed_wf_dbf.is_schedulable()) {
		//	++count_proposed_wf_dbf;
		//}
		//if (proposed_bf_den.is_schedulable()) {
		//	++count_proposed_bf_den;
		//}
		if (proposed_bf_dbf.is_schedulable()) {
			++count_proposed_bf_dbf;
		}
		*/
		if (stretch_gedf.is_schedulable()) {
			++count_stretching_gedf;
		}
		if (stretch_pedf.is_schedulable()) {
			++count_stretching_pedf;
		}
		if (decomp.is_schedulable()) {
			++count_decomp_algo;
		}
	}

	/*
	cout << "Reservation-Based [BF-RMIN-DBF] ==> Acceptance Ratio: " << count_reser_bf_dbf << "/" << num_tsets << endl;
	//	cout << "Reservation-Based [BF-RMIN-DEN] ==> Acceptance Ratio: " << count_reser_bf_den << "/" << num_tsets << endl;
	//	cout << "Reservation-Based [BF-REQUAL] ==> Acceptance Ratio: " << count_reser_requal_bf << "/" << num_tsets << endl;
	cout << "Semi-Federated ==> Acceptance Ratio: " << count_semi << "/" << num_tsets << endl;
	cout << "Baruah Algo [BF-DBF] ==> Acceptance Ratio: " << count_baruah_bf_dbf << "/" << num_tsets << endl;
	cout << "Proposed Algo [WF-DEN]  ==> Acceptance Ratio: " << count_proposed_wf_den << "/" << num_tsets << endl;
	//	cout << "Proposed Algo [WF-DBF]  ==> Acceptance Ratio: " << count_proposed_wf_dbf << "/" << num_tsets << endl;
	//	cout << "Proposed Algo [BF-DEN]  ==> Acceptance Ratio: " << count_proposed_bf_den << "/" << num_tsets << endl;
	cout << "Proposed Algo [BF-DBF]  ==> Acceptance Ratio: " << count_proposed_bf_dbf << "/" << num_tsets << endl;
	*/
	cout << "Stretching Algo [G-EDF]  ==> Acceptance Ratio: " << count_stretching_gedf << "/" << num_tsets << endl;
	cout << "Stretching Algo [P-EDF]  ==> Acceptance Ratio: " << count_stretching_pedf << "/" << num_tsets << endl;
	cout << "Decomposition Algo  ==> Acceptance Ratio: " << count_decomp_algo << "/" << num_tsets << endl;

	/*
	// Write to a result file.
	fstream out_fs(string(argv[11]) + "/result_p=" + to_string(p) + ".txt", fstream::out);
	out_fs << count_reser_bf_dbf << "/" << num_tsets << endl;
	//out_fs << count_reser_bf_den << "/" << num_tsets << endl;
	//out_fs << count_reser_requal_bf << "/" << num_tsets << endl;
	out_fs << count_semi << "/" << num_tsets << endl;
	out_fs << count_baruah_bf_dbf << "/" << num_tsets << endl;
	out_fs << count_proposed_wf_den << "/" << num_tsets << endl;
	//out_fs << count_proposed_wf_dbf << "/" << num_tsets << endl;
	//out_fs << count_proposed_bf_den << "/" << num_tsets << endl;
	out_fs << count_proposed_bf_dbf << "/" << num_tsets << endl;
	out_fs.close();
	*/
	fstream out_fs(string(argv[11]) + "/result_stretch_decomp.txt", fstream::out);
	out_fs << count_stretching_gedf << "/" << num_tsets << endl;
	out_fs << count_stretching_pedf << "/" << num_tsets << endl;
	out_fs << count_decomp_algo << "/" << num_tsets << endl;
	out_fs.close();
	return 0;
}


// MAIN FUNCTION FOR HEAVY TASK SET EXPERIMENTS.
int expr_vary_util_heavyset(int argc, char *argv[]) {
	if (argc != 8) {
		cout << "Usage: Program <NumTasksets> <MinNodes> <MaxNodes> <MinWcet> <MaxWcet> <EdgeProb> <FolderPath> !!" << endl;
		return -1;
	}

	uint_t no_tsets = stoi(string(argv[1]));
	uint_t min_size = stoi(string(argv[2]));
	uint_t max_size = stoi(string(argv[3]));
	uint_t min_wcet = stoi(string(argv[4]));
	uint_t max_wcet = stoi(string(argv[5]));
	double p = stod(string(argv[6]));
	// Path to folder that contains task sets.
	string path = string(argv[7]);

	for (int i = 1; i <= no_tsets; ++i) {
		//cout << "=============> Taskset " << i << endl;
		string tset_path = path + "/taskset" + to_string(i);
		string util_path = tset_path + "/utils.txt";
		vector<double> utils;
		fstream util_f(util_path, fstream::in);
		string util;
		while (getline(util_f, util)) {
			utils.push_back(stod(util));
		}
		util_f.close();

		for (int j = 1; j <= utils.size(); ++j) {
			DagTask task(min_size, max_size, min_wcet, max_wcet, p, utils[j-1]);
			ProposedScheduler sched(task);
			BaruahGreedyScheduler baruah(task);
			NaiveGreedyScheduler naive(task);

			string dag_file = tset_path + "/DAG_" + to_string(j) + ".txt";
			string proposed_file = tset_path + "/proposed_sched_" + to_string(j) + ".txt";
			string baruah_file = tset_path + "/baruah_sched_" + to_string(j) + ".txt";
			string naive_file = tset_path + "/naive_sched_" + to_string(j) + ".txt";
			sched.write_dagtask(dag_file);
			sched.write_schedule(proposed_file);
			baruah.write_schedule(baruah_file);
			naive.write_schedule(naive_file);
		}
	}
	
	return 0;
}


// MAIN FUNCTION FOR SINGLE HEAVY-TASK-ONLY EXPERIMENTS.
// Generate heavy dag tasks with a given edge probability and compute their
// numbers of processors needed by the proposed, baruah's and naive algorithms.
// Write the numbers of processors returned to a file.
// Write DAG files, schedule files for all generated tasks for which the proposed
// algorithm needs more cores than baruah's algorithm.
int expr_individual_heavy(int argc, char *argv[]) {
	if (argc != 8) {
		cout << "Usage: Program <NumTasks> <MinNodes> <MaxNodes> <MinWcet> <MaxWcet> <EdgeProb> <FolderPath> !!" << endl;
		return -1;
	}

	uint_t num_tasks = stoi(string(argv[1]));
	uint_t min_size = stoi(string(argv[2]));
	uint_t max_size = stoi(string(argv[3]));
	uint_t min_wcet = stoi(string(argv[4]));
	uint_t max_wcet = stoi(string(argv[5]));
	double p = stod(string(argv[6]));
	string folder = string(argv[7]);
	assert(min_size <= max_size);
	assert(min_wcet <= max_wcet);

	vector<uint_t> proposed_data(num_tasks);
	vector<uint_t> baruah_data(num_tasks);
	vector<uint_t> naive_data(num_tasks);
	
	// The number of tasks for which the proposed algorithm needs more cores.
	uint_t count_proposed_needs_more = 0;
	uint_t count_proposed_needs_less = 0;
	for (int i = 1; i <= num_tasks; ++i) {
		//cout << "Generating task " << i << "..." << endl;
		uint_t size = Common::uniform_int_gen(min_size, max_size);
		DagTask task(size, min_wcet, max_wcet, p);
		
		// Use the proposed scheduling method.
		ProposedScheduler sched(task);
		uint_t no_cores_proposed = sched.get_num_procs();
		proposed_data[i-1] = no_cores_proposed;
		
		// Use the original federated scheduling.
		NaiveGreedyScheduler naive(task);
		uint_t no_cores_naive = naive.get_schedule().size();
		naive_data[i-1] = no_cores_naive;
		
		// Use Baruah's federated scheduling.
		BaruahGreedyScheduler baruah(task);
		uint_t no_cores_baruah = baruah.get_schedule().size();
		baruah_data[i-1] = no_cores_baruah;
		
		if (no_cores_baruah < no_cores_proposed) {
			++count_proposed_needs_more;
			// Write these DAG tasks and their schedules to files.
			string dag_file = folder + "/DAG_" + to_string(i) + ".txt";
			string proposed_file = folder + "/proposed_sched_" + to_string(i) + ".txt";
			string baruah_file = folder + "/baruah_sched_" + to_string(i) + ".txt";
			string naive_file = folder + "/naive_sched_" + to_string(i) + ".txt";
			sched.write_dagtask(dag_file);
			sched.write_schedule(proposed_file);
			baruah.write_schedule(baruah_file);
			naive.write_schedule(naive_file);
			
		} else if (no_cores_baruah > no_cores_proposed) {
			++count_proposed_needs_less;
		}

		if (i % 1000 == 0) {
			cout << "Finished " << i << " tasks..." << endl;
		}
	}

	// Write the result to a file.
	fstream of(folder + "/result.txt", fstream::out);
	for (int i = 0; i < num_tasks; ++i) {
		string line = to_string(proposed_data[i]) + "\t" + to_string(baruah_data[i]) + "\t" + to_string(naive_data[i]);
		of << line << endl;
	}
	of.close();
	
	cout << "#Times the Proposed Algo needs MORE cores: " << count_proposed_needs_more << "/" << num_tasks << endl;
	cout << "#Times the Proposed Algo needs LESS cores: " << count_proposed_needs_less << "/" << num_tasks << endl;
	
	return 0;
}

int test_extra_preemptions(int argc, char *argv[]) {
	/*
	matrix_t adj_matrix{{0, 0, 0, 1, 0, 0, 0, 1, 0, 1},
						{0, 0, 0, 0, 0, 1, 0, 1, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
	vector<uint_t> wcets{18, 9, 9, 9, 9, 18, 18, 18, 7, 7};
	uint_t period = 44, deadline = 44;
	DagTask task(adj_matrix, wcets, period, deadline);
	*/
	if (argc != 2) {
		cout << "Usage: Program <NumSubtasks> !!" << endl;
		return -1;
	}
	uint_t size = stoi(argv[1]);
	
	DagTask task(size, 50, 100, 0.2);
	ProposedScheduler sched(task);
	auto subtasks = sched.get_schedule();
	uint_t no_preempts = 0;
	for (const Subtask& s : subtasks) {
		const list<Fragment>& frags = s.frags;
		uint_t prev_finish = frags.front().start_time;
		for (const Fragment f : frags) {
			if (f.start_time > prev_finish) {
				++no_preempts;
			}
			prev_finish = f.start_time + f.exe_time;
		}
	}

	cout << "Number of extra preemptions: " << no_preempts << endl;
	
	// Write the DAG task and the schedules to files.
	string dag_file = "tmp/DAG.txt";
	string proposed_file = "tmp/proposed_sched.txt";
	sched.write_dagtask(dag_file);
	sched.write_schedule(proposed_file);
	
	return 0;
}

// Read task utilizations for a task set from a file.
// Generate DAG tasks with those utilizations and run scheduling algorithms.
int main9(int argc, char *argv[]) {
	if (argc != 5) {
		cout << "Usage: Program <UtilsFilePath> <NumSubtasks> <EdgeProb> <OutPath> !!" << endl;
		return -1;
	}

	// Read utilizations.
	vector<double> utils;
	fstream util_f(argv[1], fstream::in);
	string util;
	while (getline(util_f, util)) {
		utils.push_back(stod(util));
	}
	util_f.close();

	// Generating DAG tasks with the read utilizations.
	// Run scheduling algorithms for each task.
	uint_t size = stoi(string(argv[2]));
	uint_t wmin = 1, wmax = 20;
	double p = stod(argv[3]);
	string out_folder(argv[4]); // E.g.: ./tmp/taskset1/
	for (int i = 1; i <= utils.size(); ++i) {
		cout << "======== Generating task " << i << " ..." << endl;
		DagTask task(size, wmin, wmax, p, utils[i-1]);
		ProposedScheduler sched(task);
		BaruahGreedyScheduler baruah(task);
		NaiveGreedyScheduler naive(task);
		
		// Write the DAG task and the schedules to files.
		string dag_file = out_folder + "/DAG_" + to_string(i) + ".txt";
		string proposed_file = out_folder + "/proposed_sched_" + to_string(i) + ".txt";
		string baruah_file = out_folder + "/baruah_sched_" + to_string(i) + ".txt";
		string naive_file = out_folder + "/naive_sched_" + to_string(i) + ".txt";
		sched.write_dagtask(dag_file);
		sched.write_schedule(proposed_file);
		baruah.write_schedule(baruah_file);
		naive.write_schedule(naive_file);
	}
	
	return 0;
}

// Read a DAG file and run scheduling algorithms for that DAG.
// This is for debugging a couple of tasks for which the proposed algorithm
// needs more processors. Fixed on Sept. 27, 2019.
int main8(int argc, char *argv[]) {
	if (argc != 2) {
		cout << "Usage: Program <DagFile> !!" << endl;
		return -1;
	}

	string file_path(argv[1]);
	vector<uint_t> wcets;
	matrix_t adj_matrix;
	uint_t period, deadline;
	read_dag_from_file(file_path, wcets, adj_matrix, period, deadline);
	DagTask task(adj_matrix, wcets, period, deadline);
	ProposedScheduler sched(task);
	string proposed_file = "tmp/proposed_sched_133639_new.txt";
	sched.write_schedule(proposed_file);
	return 0;
}

// Run a large number of tasks and check if there are any tasks for which the proposed
// needs more processors than the improved greedy algorithm. Only write the DAGs and
// the schedules to files for such tasks.
int main7(int argc, char *argv[]) {
	if (argc != 5) {
		cout << "Usage: Program <NumTasks> <NumSubtasks> <EdgeProb> <FolderPath> !!" << endl;
		return -1;
	}

	uint_t num_tasks = stoi(string(argv[1]));
	uint_t size = stoi(string(argv[2]));
	uint_t wmin = 1, wmax = 20;
	double p = stod(string(argv[3]));
	string folder = string(argv[4]);

	// The number of tasks for which the proposed algorithm needs more cores.
	uint_t count_proposed_needs_more = 0;
	uint_t count_proposed_needs_less = 0;
	for (int i = 1; i <= num_tasks; ++i) {
		//cout << "Generating task " << i << "..." << endl;
		DagTask task(size, wmin, wmax, p);
		
		// Use the proposed scheduling method.
		ProposedScheduler sched(task);
		uint_t no_cores_proposed = sched.get_num_procs();
		
		// Use the original federated scheduling.
		NaiveGreedyScheduler naive(task);
		
		// Use Baruah's federated scheduling.
		BaruahGreedyScheduler baruah(task);
		uint_t no_cores_baruah = baruah.get_schedule().size();

		if (no_cores_baruah < no_cores_proposed) {
			++count_proposed_needs_more;
			// Write these DAG tasks and their schedules to files.
			string dag_file = folder + "/DAG_" + to_string(i) + ".txt";
			string proposed_file = folder + "/proposed_sched_" + to_string(i) + ".txt";
			string baruah_file = folder + "/baruah_sched_" + to_string(i) + ".txt";
			string naive_file = folder + "/naive_sched_" + to_string(i) + ".txt";
			sched.write_dagtask(dag_file);
			sched.write_schedule(proposed_file);
			baruah.write_schedule(baruah_file);
			naive.write_schedule(naive_file);
			
		} else if (no_cores_baruah > no_cores_proposed) {
			++count_proposed_needs_less;
		}

		if (i % 1000 == 0) {
			cout << "Finished " << i << " tasks..." << endl;
		}
	}
	cout << "#Times the Proposed Algo needs MORE cores: " << count_proposed_needs_more << "/" << num_tasks << endl;
	cout << "#Times the Proposed Algo needs LESS cores: " << count_proposed_needs_less << "/" << num_tasks << endl;
	
	return 0;
}


// Debugging the improved greedy algorithm
int main6(int argc, char *argv[]) {
	// Read task from 10nodes_p=0.3/DAG_4.txt
	matrix_t adj_matrix;
	vector<uint_t> wcets;
	uint_t period, deadline;
	string path = "data/new_algo/10nodes_p=0.3/DAG_4.txt";
	read_dag_from_file(path, wcets, adj_matrix, period, deadline);

	// Create DAG task object.
	DagTask task(adj_matrix, wcets, period, deadline);
	BaruahGreedyScheduler baruah(task);
	uint_t num_cores = baruah.get_schedule().size();
	string out_file = "debugging_baruah_sched_" + to_string(num_cores) + "cores.txt";
	baruah.write_schedule(out_file);
	return 0;
}

// Generate random DAG tasks, perform scheduling and write the resulting schedules to files.
int main5(int argc, char *argv[]) {
	if (argc != 5) {
		cout << "Usage: Program <Num_Tasks> <Num_Subtasks> <EdgeProb> <FolderPath> !!" << endl;
		return -1;
	}

	uint_t num_tasks = stoi(string(argv[1]));
	uint_t size = stoi(string(argv[2]));
	uint_t wmin = 1, wmax = 20;
	double p = stod(string(argv[3]));
	string folder = string(argv[4]);
		
	for (int i = 1; i <= num_tasks; ++i) {
		DagTask task(size, wmin, wmax, p);
		//print_dag_task(task);
		
		// Use the proposed scheduling method.
		ProposedScheduler sched(task);
		//print_proposed_schedule(sched);
		
		// Use the original federated scheduling.
		NaiveGreedyScheduler naive(task);
		//print_naive_schedule(naive, task);
		
		// Use Baruah's federated scheduling.
		BaruahGreedyScheduler baruah(task);
		//print_baruah_schedule(baruah, task);
		
		// Write the DAG task and the schedules to files.
		string dag_file = folder + "/DAG_" + to_string(i) + ".txt";
		string proposed_file = folder + "/proposed_sched_" + to_string(i) + ".txt";
		string baruah_file = folder + "/baruah_sched_" + to_string(i) + ".txt";
		string naive_file = folder + "/naive_sched_" + to_string(i) + ".txt";
		sched.write_dagtask(dag_file);
		sched.write_schedule(proposed_file);
		baruah.write_schedule(baruah_file);
		naive.write_schedule(naive_file);
	}
	
	return 0;
}

void read_dag_from_file(string dag_file, vector<uint_t> &wcets, matrix_t &adj_matrix, uint_t &period, uint_t &deadline) {
	fstream fs(dag_file, fstream::in);
	// Read number of subtasks.
	string no_subtasks_str;
	getline(fs, no_subtasks_str);
	uint_t no_subtasks = stoi(no_subtasks_str);
	
	// Read the WCETs of the subtasks.
	string wcets_str;
	getline(fs, wcets_str);
	stringstream wcets_ss(wcets_str);
	uint_t temp;
	while (wcets_ss >> temp) {
		wcets.push_back(temp);
	}
	
	// Read the task's parameters.
	string param_str;
	getline(fs, param_str);
	stringstream param_ss(param_str);
	uint_t work, span;
	param_ss >> work;
	param_ss >> span;
	param_ss >> deadline;
	param_ss >> period;
	
	// Read the critical-path.
	string span_str;
	getline(fs, span_str);
	stringstream span_ss(span_str);
	vector<uint_t> crit_path;
	while (span_ss >> temp) {
		crit_path.push_back(temp);
	}
	
	// Read adjacency list of the DAG.
	list_t adj_list(no_subtasks); // vector of lists.
	string line;
	while (getline(fs, line)) {
		stringstream adj_ss(line);
		uint_t sub_id;
		adj_ss >> sub_id;
		while (adj_ss >> temp) {
			adj_list[sub_id].push_back(temp);
		}
	}
	
	fs.close();
	
	// Form the adjacency matrix from the adjacency list.
	adj_matrix.resize(no_subtasks);
	for (uint_t i = 0; i < no_subtasks; ++i) {
		adj_matrix[i].resize(no_subtasks, 0);
		for (uint_t j : adj_list[i]) {
			adj_matrix[i][j] = 1;
		}
	}	
}

// Read DAGs from files and run the scheduling algorithms.
int main4(int argc, char *argv[]) {
	if (argc != 3) {
		cout << "Usage: Program <Folder> <NumDags> !!" << endl;
		return -1;
	}

	uint_t no_dags = stoi(argv[2]);

	for (uint_t i = 1; i <= no_dags; ++i) {
		string dag_file = string(argv[1]) + "/DAG_" + to_string(i) + ".txt";
		matrix_t adj_matrix;
		vector<uint_t> wcets;
		uint_t period, deadline;
		read_dag_from_file(dag_file, wcets, adj_matrix, period, deadline);
		
		// Create the task instance and call the schedulers.
		DagTask task(adj_matrix, wcets, period, deadline);
		ProposedScheduler sched(task);
		BaruahGreedyScheduler baruah(task);
		NaiveGreedyScheduler naive(task);
		
		string proposed_file = string(argv[1]) + "/proposed_sched_" + to_string(i) + ".txt";
		string baruah_file = string(argv[1]) + "/baruah_sched_" + to_string(i) + ".txt";
		string naive_file = string(argv[1]) + "/naive_sched_" + to_string(i) + ".txt";
		sched.write_schedule(proposed_file);
		baruah.write_schedule(baruah_file);
		naive.write_schedule(naive_file);
	}
	
	return 0;
}


void print_proposed_schedule(const ProposedScheduler& sched) {
	cout << "========= Proposed Schedule ===========" << endl;
	//	const vector<Subtask>& subtasks = sched.get_schedule();
	auto subtasks = sched.get_schedule();
	for (const Subtask& sub : subtasks) {
		for (list<Fragment>::const_iterator it = sub.frags.cbegin(); it != sub.frags.cend(); ++it) {
			cout << "Subtask " << it->sub_id << ". Start:  " << it->start_time << ". Exe time: " <<
				it->exe_time << ". Deadline: " << it->deadline << ". Proc: " << it->proc_id << "." << endl;
		}
	}
}

void print_dag_task(const DagTask& task) {
	cout << "========= DAG Structure ==========" << endl;
	cout << "Node lengths: ";
	for (uint_t wcet : task.get_node_lengths()) {
		cout << wcet << ", ";
	}
	auto crit_path = task.get_crit_path();
	cout << "Critical-path: ";
	for (uint_t sub_id : crit_path) {
		cout << sub_id << " ";
	}

	cout << ". Work: " << task.get_work() << ". Span: " << task.get_span() << ". Deadline: " << task.get_deadline() << endl;
	cout << "DAG: ";
	for (uint_t i = 0; i < task.get_adj_list().size(); ++i) {
		cout << "Subtask " << i << ": ";
		for (uint_t j : task.get_adj_list()[i]) {
			cout << j << ", ";
		}
		cout << endl;
	}
}

void print_naive_schedule(const NaiveGreedyScheduler& naive, const DagTask& task) {
	cout << "========== Naive Greedy Schedule ===========" << endl;
	const vector<Processor>& procs = naive.get_schedule();
	
	for (const Processor& proc : procs) {
		cout << "Processor " << proc.id << " ==> ";
		for (FinishedSubtask t : proc.subtasks) {
			cout << "Subtask: " << t.id << ", start: " << t.start_time << ", wcet: "
				 << task.get_node_lengths()[t.id] << ". ";
		}
		cout << endl;
	}
}

void print_baruah_schedule(const BaruahGreedyScheduler& baruah, const DagTask& task) {
	cout << "========== Baruah Greedy Schedule ===========" << endl;
	const vector<Processor>& baruah_procs = baruah.get_schedule();

	for (const Processor& proc : baruah_procs) {
		cout << "Processor " << proc.id << " ==> ";
		for (FinishedSubtask t : proc.subtasks) {
			cout << "Subtask: " << t.id << ", start: " << t.start_time << ", wcet: "
				 << task.get_node_lengths()[t.id] << ". ";
		}
		cout << endl;
	}
}


int main3() {
	/*
	matrix_t adj_matrix{{0,1,1,1,1,1,1,1,1,1},
						{0,0,1,1,1,1,1,1,1,1},
						{0,0,0,1,1,1,1,0,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,0}};
	vector<uint_t> wcets{3, 9, 9, 9, 9, 9, 9, 16, 16, 16};
	uint_t period = 32, deadline = 32;
	
	// Example DAG task in the submission to RTAS 2020.
	matrix_t adj_matrix{{0,0,0,1,0,0,0,1,0,1},
						{0,0,0,0,0,1,0,1,0,0},
						{0,0,0,0,0,0,0,1,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,1,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,1,0,0},
						{0,0,0,0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0,0,0,1},
						{0,0,0,0,0,0,0,0,0,0}};
	vector<uint_t> wcets{18, 9, 9, 9, 9, 18, 18, 18, 7, 7};
	uint_t period = 44, deadline = 44;
	
	matrix_t adj_matrix{{0,1,1,0,0,0,0},
						{0,0,0,0,1,0,0},
						{0,0,0,0,1,1,1},
						{0,0,0,0,1,0,1},
						{0,0,0,0,0,0,1},
						{0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0}};
	vector<uint_t> wcets{2, 5, 3, 2, 2, 4, 3};
	uint_t period = 13, deadline = 13;
	
	// Example from the Qamhieh's stretching paper.
	matrix_t adj_matrix{{0,0,0,1,0,0,0},
						{0,0,0,1,0,0,0},
						{0,0,0,0,0,1,0},
						{0,0,0,0,0,1,1},
						{0,0,0,0,0,0,1},
						{0,0,0,0,0,0,0},
						{0,0,0,0,0,0,0}};
	vector<uint_t> wcets{3, 3, 2, 1, 2, 2, 1};
	uint_t period = 10, deadline = 10;
	*/	
	// Example from the Jiang's decomposition paper.
	matrix_t adj_matrix{{0,1,1,1,0,0,0,0},
						{0,0,0,0,0,0,0,1},
						{0,0,0,0,1,0,1,0},
						{0,0,0,0,0,0,1,0},
						{0,0,0,0,0,1,1,0},
						{0,0,0,0,0,0,0,1},
						{0,0,0,0,0,0,0,1},
						{0,0,0,0,0,0,0,0}};
	vector<uint_t> wcets{1, 3, 4, 8, 5, 4, 4, 1};
	uint_t period = 25, deadline = 25;
	
	DagTask task(adj_matrix, wcets, period, deadline);

	ProposedScheduler sched(task);
	//print_proposed_schedule(sched);

	NaiveGreedyScheduler naive(task);
	//print_naive_schedule(naive, task);
	
	// Use Baruah's federated scheduling.
	BaruahGreedyScheduler baruah(task);
	//print_baruah_schedule(baruah, task);
	
	// Write the DAG task and the schedules to files.
	//	string dag_file = "./data/DAG_failed.txt";
	//	string proposed_file = "./data/proposed_sched_failed.txt";
	//	string baruah_file = "./data/baruah_sched_failed.txt";
	//	string naive_file = "./data/naive_sched_failed.txt";
	string dag_file = "./tmp/DAG_example.txt";
	string proposed_file = "./tmp/proposed_sched_exp.txt";
	string baruah_file = "./tmp/baruah_sched_exp.txt";
	sched.write_dagtask(dag_file);
	sched.write_schedule(proposed_file);
	baruah.write_schedule(baruah_file);
	//	naive.write_schedule(naive_file);

	Taskset tset(2);
	tset.add_heavy_task(task);
	//StretchingAlgo stretch(tset);
	//stretch.is_schedulable();
	DecompAlgo dcom(tset);
	dcom.is_schedulable();
	
	return 0;
}

int main1() {
	DagGen dag;
	dag.test();
	return 0;
}

int main2(int argc, char *argv[]) {
	if (argc != 3) {
		cout << "Usage: Program Size Prob" << endl;
		return -1;
	}
	unsigned int size = stoi(argv[1]);	
	double p = stod(argv[2]);;

	// Generate adjacency matrix.
	DagGen dag;
	vector<vector<unsigned int> > am;
	dag.gen_erdos_matrix(am, size, p);

	// Generate subtasks' lengths.
	vector<unsigned int> wcets;
	unsigned int work = dag.gen_node_lengths(wcets, size, 1, 20);

	// Compute critical path length.
	unsigned int span_old = dag.calc_span_old(am, wcets);

	// Test new span computation.
	vector<pair<unsigned int, int>> dist;
	vector<unsigned int> crit_path;
	unsigned int span_new = dag.calc_span(am, wcets, dist, crit_path);

	cout << "Work: " << work << ". Old span: " << span_old << ". New span: " << span_new << endl;
	for (unsigned int node : crit_path) {
		cout << "(" << node << ", " << wcets[node] << ") ";
	}
	cout << endl;
	
	cout << "Nodes lengths: ";
	for (unsigned long a : wcets) {
		cout << a << " ";
	}
	cout << endl;

	for (int i = 0; i < size; ++i) {
		unsigned int temp = dag.calc_subdag_work(i, am, wcets);
		cout << "Total work node " << i << ": " << temp << endl;
	}
	
	for (vector<unsigned>& v : am) {
		for (unsigned a : v) {
			cout << a << " ";
		}
		cout << endl;
	}

	return 0;
}
