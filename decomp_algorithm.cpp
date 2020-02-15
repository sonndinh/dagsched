#include <queue>
#include <algorithm>
#include "decomp_algorithm.h"

DecompAlgo::DecompAlgo(const Taskset& taskset) : tset_{taskset} {}

// Apply Theorem 3 in Jiang et al.'s paper.
bool DecompAlgo::is_schedulable() {
	for (const DagTask& t : tset_.get_heavy_tasks()) {
		decompose(&t);
	}

	for (const SequentialTask& t : tset_.get_light_tasks()) {
		decompose(&t);
	}

	// Total utilization of the task set.
	double total_util = 0;
	// Maximum ratio of L_i/T_i of all tasks.
	double gamma = 0.0;

	for (const DagTask& t : tset_.get_heavy_tasks()) {
		total_util += t.get_utilization();
		gamma = max(gamma, (double)t.get_span()/t.get_period());
	}
	
	for (const SequentialTask& t : tset_.get_light_tasks()) {
		total_util += t.get_utilization();
		gamma = max(gamma, (double)t.get_span()/t.get_period());
	}

	// Maximum structure characteristic values of all tasks.
	double omega = 0.0;
	for (const TimingDiagram& d : diagrams_) {
		omega = max(omega, d.omega);
	}
	
	// Test the condition in Theorem 3.
	double tmp = (total_util - gamma)/(1/omega - gamma);
	if (tset_.get_no_procs() >= tmp)
		return true;
	return false;
}

// Implement the segmentation step to decompose a given DAG task
// (heavy or light) into a set of sequential tasks.
void DecompAlgo::decompose(const Task *task) {
	uint_t work = task->get_work();
	uint_t span = task->get_span();
	uint_t idx = compute_timing_diagram(task);
	TimingDiagram& diagram = diagrams_[idx];
	const vector<uint_t>& wcets = task->get_node_lengths();
	priority_queue<NodeWindow, vector<NodeWindow>, IncreasingFinish> queue(diagram.subtasks.begin(), diagram.subtasks.end());
	
	// Phase 1: assign vertices that only cover a single segment.
	vector<NodeWindow> tmp;
	while (!queue.empty()) {
		const NodeWindow& node = queue.top();
		bool is_assigned = false;
		for (DcomSegment& seg : diagram.segments) {
			if (seg.start == node.ready && seg.end == node.finish) {
				// Add this node to this segment.
				seg.id_to_len[node.id] = node.len;
				seg.work += node.len;
				queue.pop();
				is_assigned = true;
				break;
			}
		}
		if (is_assigned == false) {
			tmp.push_back(node);
			queue.pop(); // Temporarily remove from the queue.
		}
	}
	// Insert the unassigned nodes back into the queue.
	for (NodeWindow& n : tmp) {
		queue.push(n);
	}
	

	// Phase 2: assign the remaining nodes to light segments.
	// Implemented exactly as the pseudocode in the paper.
	for (DcomSegment& seg : diagram.segments) {
		// Ignore heavy segments.
		if ((double)seg.work/(seg.end - seg.start) >= (double)work/span)
			continue;

		tmp.clear();
		while (!queue.empty()) {
			const NodeWindow& node = queue.top();
			if (!(node.finish <= seg.start || node.ready >= seg.end)) {
				double lhs = (double)(node.len + seg.work)/(seg.end - seg.start);
				if (lhs <= (double)work/span) {
					// Add the entire node to the segment.
					seg.id_to_len[node.id] = node.len;
					seg.work += node.len;
					queue.pop();
				} else { // Split the node.
					double insert_len = ((double)work/span)*(seg.end - seg.start) - (double)seg.work;
					if (insert_len < 0) {
						cout << "ERROR: Inserting node " << node.id << " with negative length " << insert_len << " !!" << endl;
					}
					seg.id_to_len[node.id] = insert_len;
					seg.work += insert_len;
					NodeWindow remain = node;
					remain.len -= insert_len;
					//remain.ready += insert_len;
					queue.pop();
					queue.push(remain);
					break;
				}
			} else { // Temporarily remove the node from queue.
				tmp.push_back(node);
				queue.pop();
			}
		}

		for (NodeWindow& node : tmp) {
			queue.push(node);
		}
	}

	/*	
	// Phase 2: assign the remaining nodes to light segments.	
	// NOTE: this implementation of phase 2 does not strictly follow Algorithm 1 in Jiang et al.'s paper.
	// Specifically, a node that can be inserted entirely into a light segment without making the segment
	// becoming heavy, and has length greater than the length of the targeted segment is split.
	// Only a part of it with length equal to the segment length is inserted.
	// Phases 1 and 3 are kept the same.
	for (DcomSegment& seg : diagram.segments) {
		// Ignore heavy segments in this phase.
		if ((double)seg.work/(seg.end - seg.start) >= (double)work/span)
			continue;

		tmp.clear();
		while (!queue.empty()) {
			const NodeWindow& node = queue.top();
			if (!(node.finish <= seg.start || node.ready >= seg.end)) {
				// If the node covers the segment, add it to the segment.
				double lhs = (double)(node.len + seg.work)/(seg.end - seg.start);
				if (lhs <= (double)work/span) {
					// The segment is light even if the entire node is added.
					if (node.len <= seg.end - seg.start) {
						// Add the entire node to the segment.
						//cout << "1 -- Inserting node " << node.id << " with length " << node.len << endl;
						seg.id_to_len[node.id] = node.len;
						seg.work += node.len;
						queue.pop();
					} else {
						// Split the node and add a part with length
						// equal to the segment length.
						//cout << "2 -- Inserting node " << node.id << " with length " << seg.end - seg.start << endl;
						seg.id_to_len[node.id] = seg.end - seg.start;
						seg.work += seg.end - seg.start;
						// Temporarily save the remain for inserting back to
						// the queue later.
						NodeWindow remain = node;
						remain.len -= seg.end - seg.start;
						//remain.ready += seg.end - seg.start;
						queue.pop();
						tmp.push_back(remain);
					}
				} else { // Have to split the node.
					double max_insert = ((double)work/span)*(seg.end - seg.start) - (double)seg.work;
					// Can only add a part with length of at most the
					// length of the segment.
					max_insert = min(max_insert, (double)(seg.end - seg.start));
					if (max_insert < 0) {
						cout << "ERROR: Node " << node.id << " is split and inserted " << max_insert << " unit!!" << endl;
					}
					//cout << "3 -- Inserting node " << node.id << " with length " << max_insert << endl;
					seg.id_to_len[node.id] = max_insert;
					seg.work += max_insert;
					NodeWindow remain = node;
					remain.len -= max_insert;
					//remain.ready += max_insert;
					queue.pop();
					if (remain.len > 0) {
						tmp.push_back(remain);
					}
				}

				// Move to the next segment if this segment is done.
				if ((double)seg.work/(seg.end - seg.start) >= (double)work/span)
					break;
			} else { // This node does not cover this segment.
				tmp.push_back(node);
				queue.pop();
			}
		}
		
		for (NodeWindow& n : tmp) {
			//cout << "(node " << n.id << ", len:" << n.len << ") ";
			queue.push(n);
		}
		//		cout << endl;
	}
	*/
	
	// Phase 3: assign the remaining nodes arbitrarily to the segments
	// covered by each node.
	while (!queue.empty()) {
		NodeWindow node = queue.top();
		queue.pop();
		for (DcomSegment& seg : diagram.segments) {
			if (!(node.finish <= seg.start || node.ready >= seg.end)) {
				if ((seg.id_to_len.find(node.id) != seg.id_to_len.end() &&
					 seg.id_to_len[node.id] < seg.end - seg.start)) {
					uint_t added_len = min(node.len, seg.end - seg.start - seg.id_to_len[node.id]);
					seg.id_to_len[node.id] += added_len;
					seg.work += added_len;
					node.len -= added_len;
					//node.ready += added_len;
					if (node.len == 0) {
						break;
					}
				} else if (seg.id_to_len.find(node.id) == seg.id_to_len.end()) {
					uint_t added_len = min(node.len, seg.end - seg.start);
					seg.id_to_len[node.id] = added_len;
					seg.work += added_len;
					node.len -= added_len;
					//node.ready += added_len;
					if (node.len == 0) {
						break;
					}
				}
			}
		}

		// Debugging.
		if (node.len > 0) {
			cout << "ERROR: Node " << node.id << " was not assigned successfully in Phase 3 !!" << endl;
		}
	}

	// Compute Omega for this task.
	uint_t heavy_work = 0, light_len = 0;
	for (const DcomSegment& seg : diagram.segments) {
		if ((double)seg.work/(seg.end - seg.start) > (double)work/span) {
			heavy_work += seg.work;
		} else {
			light_len += seg.end - seg.start;
		}
	}
	diagram.omega = (double)heavy_work/work + (double)light_len/span;

	// For testing
	//	cout << " ====== Decomposed task: " << endl;
	//	print_segments(idx);
}

// Compute the timing diagram for a given task that contains
// the earliest ready time and latest finish time of each node. 
// Return the index of the diagram for this task.
uint_t DecompAlgo::compute_timing_diagram(const Task* task) {
	TimingDiagram diagram;
	uint_t size = task->get_size();
	diagram.subtasks.resize(size);
	const vector<uint_t>& wcets = task->get_node_lengths();
	const matrix_t& adj_matrix = task->get_adj_matrix();
	list<uint_t> topo = topological_sort(adj_matrix);
	vector<uint_t> sorted_nodes(topo.begin(), topo.end());
	unordered_set<uint_t> ready_times;

	// Compute the earliest ready times for the subtasks.
	// Traverse the nodes in their topological order.
	for (uint_t k = 0; k < size; ++k) {
		uint_t node_id = sorted_nodes[k];
		bool is_source = true;
		for (uint_t i = 0; i < size; ++i) {
			if (adj_matrix[i][node_id] == 1) {
				is_source = false;
				break;
			}
		}

		NodeWindow tmp;
		tmp.id = node_id;
		tmp.len = wcets[node_id];
		if (is_source) {
			// Source nodes are ready at time 0.
			tmp.ready = 0;
		} else {
			uint_t ready = 0;
			// Traverse through all predecessors of this node.
			for (uint_t j = 0; j < k; ++j) {
				uint_t pred_id = sorted_nodes[j];
				if (adj_matrix[pred_id][node_id] == 1) {
					// The @ready of node pred_id has already been computed.
					ready = max(ready, wcets[pred_id] + diagram.subtasks[pred_id].ready);
				}
			}
			tmp.ready = ready;
		}
		ready_times.insert(tmp.ready);
		diagram.subtasks[node_id] = tmp;
	}

	// Compute the latest finish times for all subtasks.
	for (uint_t k = 0; k < size; ++k) {
		uint_t node_id = sorted_nodes[k];
		// Traverse through all successors of this node.
		uint_t finish = task->get_period();
		for (uint_t j = k + 1; j < size; ++j) {
			uint_t succ_id = sorted_nodes[j];
			if (adj_matrix[node_id][succ_id] == 1) {
				finish = min(finish, diagram.subtasks[succ_id].ready);
			}
		}
		if (k == size - 1) {
			diagram.subtasks[node_id].finish = task->get_span();
		} else {
			diagram.subtasks[node_id].finish = finish;
		}
	}

	vector<uint_t> segment_times(ready_times.begin(), ready_times.end());
	sort(segment_times.begin(), segment_times.end());
	diagram.segments.resize(segment_times.size());
	
	for (uint_t i = 0; i < segment_times.size(); ++i) {
		diagram.segments[i].start = segment_times[i];
		if (i == segment_times.size() - 1) {
			diagram.segments[i].end = task->get_span();
		} else {
			diagram.segments[i].end = segment_times[i+1];
		}
	}
		
	diagrams_.push_back(diagram);
	
	// For testing
	//	cout << "======== Timing diagram: " << endl;
	//	print_timing_diagram(diagrams_.size() - 1);
	
	return (diagrams_.size() - 1);
}

list<uint_t> DecompAlgo::topological_sort(const matrix_t& adj) {
	list<uint_t> ret;
	uint_t size = adj.size();
	unordered_set<uint_t> visited;
	
	for (uint_t i = 0; i < size; ++i) {
		if (visited.find(i) == visited.end()) {
			dfs_visit(adj, i, ret, visited);
		}
	}

	return ret;
}


void DecompAlgo::dfs_visit(const matrix_t& adj, uint_t i, list<uint_t>& ret, unordered_set<uint_t>& visited) {
	visited.insert(i);
	for (uint_t j = 0; j < adj.size(); ++j) {
		if (adj[i][j] == 1 && visited.find(j) == visited.end()) {
			dfs_visit(adj, j, ret, visited);
		}
	}
	ret.push_front(i);
}

// Print the timing diagram of a given task.
void DecompAlgo::print_timing_diagram(uint_t idx) {
	const TimingDiagram& diagram = diagrams_[idx];
	uint_t size = diagram.subtasks.size();
	
	for (uint_t i = 0; i < size; ++i) {
		cout << "Node " << diagram.subtasks[i].id << " -- Ready: "
			 << diagram.subtasks[i].ready << ". Finish: "
			 << diagram.subtasks[i].finish << endl;
	}
}

void DecompAlgo::print_segments(uint_t idx) {
	const TimingDiagram& diagram = diagrams_[idx];
	uint_t no_segments = diagram.segments.size();

	for (uint_t i = 0; i < no_segments; ++i) {
		cout << "Segment " << i << " -- Start: " << diagram.segments[i].start
			 << ". End: " << diagram.segments[i].end
			 << ". Total work: " << diagram.segments[i].work
			 << " ==> ";
		for (auto p : diagram.segments[i].id_to_len) {
			cout << "(id:" << p.first << ", len:" << p.second << ") ";
		}
		cout << endl;
	}
}
