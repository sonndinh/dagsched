// Generate a random connected directed-acyclic-graph.
// Original author: David Ferry.
// Modified by: Son Dinh. 


#include <iostream>
#include <list>
#include <cstdlib>
#include <vector>
#include <stack>
#include <unordered_set>
#include <queue>
#include <numeric>
#include "rand_dag.h"

using namespace std;

//Suppose a graph has no cycles in it, and we want to check if adding a
//particular edge introduces a cycle. Since there are no previous cycles
//in the graph this new edge must be a part of any introduced cycle.
//Thus, if we want to add a new edge originating in "start" and 
//terminating in "end" then we need only do a DFS from "end" to see
//if "start" is reachable. 

//If this method returns "false" then a cycle exists, if it returns
//"true" then no cycle exists
bool DagGen::acyclic_check(matrix_t &am, uint_t size, uint_t start, uint_t end) {
	if( start == end){
		return false;
	}
	
	vector<bool> visited (size,false);
	uint_t current;
	
	//DFS code from here to end of routine
	list<uint_t> stack;
	stack.push_front(end);
	
	while(!stack.empty()) {
		current = stack.front();
		stack.pop_front();
		visited[current] = true;
	
		//If the start node is reachable from any node on the stack then there
		//is a cycle in the graph
		if(am[current][start] == 1) {
			return false;
		}

		//If not we push all unvisited children onto the stack
		for(uint i = 0; i < size; i++) {
			if(am[current][i] == 1 && visited[i] == false) {
				stack.push_front(i);
			}
		}
	} 

	//If we've searched the whole graph and not gotten back to start then
	//there is no cycle in the graph. 
	return true;
}

//A digraph is weakly connected if a DFS along the undirected analog can reach
//every node in the graph
bool DagGen::weak_conn_test(matrix_t &am) {
	uint_t size = am.size();
	vector<bool> visited (size, false);
	list<uint_t> stack;
	stack.push_front(0);

	uint_t current;

	//DFS traversal with undirected edges
	while(!stack.empty()) {
		current = stack.front();
		stack.pop_front();
		visited[current] = true;
		
		//We push all unvisited children onto the stack
		for(uint_t i = 0; i < size; i++) {
			//Forward edges
			if(am[current][i] == 1 && visited[i] == false) {
				stack.push_front(i);
			}
			//Backward edges
			if(am[i][current] == 1 && visited[i] == false) {
				stack.push_front(i);
			}
		}
	}
	
	//Test to see if every node was visited or not	
	bool connected = true;
	for(vector<bool>::iterator it = visited.begin(); it != visited.end(); it++) {
		if( *it == false ) {
			connected = false;
		}
	}
	
	return connected;
}

//Generates an adjacency matrix representation of a DAG
void DagGen::gen_adj_matrix(matrix_t &am, uint_t size) {

	//Get size of DAG to generate
	uint_t n = size;

	//Declare an adjacency matrix and initialize it to zeroes
	am.resize(n);
	for(uint_t i = 0; i < n; i++) {
		am[i].resize(n);
	}

	for(uint_t i = 0; i < n; i++) {
		for(uint_t j = 0; j < n; j++) {
			am[i][j] = 0;
		}
	}

	#ifdef VERBOSE
	cout << "Generating directed graph of size " << n << endl;
	#endif

	uint_t node;
	uint_t target;
	bool graph_is_connected = false;

	//We continue to add edges (i.e. loop) until we have a DAG that weakly
	//connects all vertices
	while(!graph_is_connected) {
		//		node = rand()%n;
		//		target = rand()%n;

		// Son -- begin
		node = Common::uniform_int_gen(0, n-1);
		target = Common::uniform_int_gen(0, n-1);
		// Son -- end
		
		//we've already looked at this edge
		if(am[node][target] != 0) {
			continue;
		}
		
		switch (am[node][target]) {
			//Case: we haven't looked at this edge yet
		case 0:
			//Check to see if adding this edge adds a cycle
			if(acyclic_check(am, n, node, target)) {
				//If no cycle is introduced, add the edge to the graph
				am[node][target] = 1;
				
				//If the graph is now weakly connected then we terminate
				if(weak_conn_test(am)) {
					graph_is_connected = true;
				}	
				
			} else {
				//If a cycle is introduced, mark it as such
				am[node][target] = 2;
			}	
			break;
	
			//Case: this edge already exists
		case 1:
			break;
			
			//Case: we already looked at this edge and discarded it
		case 2:
			break;
		default:
			cout << "Found bad value in adjacency matrix" << endl;
			abort();
		}
	}
	
	/** For debugging **
		cout << "Resultant graph:" << endl;
		for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
		cout << am[i][j] << " ";
		}
		cout << endl;
		}
	*/
	
}

//A digraph is weakly connected if a DFS along the undirected analog can reach
//every node in the graph
int DagGen::weak_conn_add(matrix_t &am){
	
	uint_t size = am.size();
	vector<bool> visited (size,false);
	list<uint_t> stack;
	stack.push_front(0);	
	uint_t current;

	//DFS traversal with undirected edges
	while (!stack.empty()) {
		current = stack.front();
		stack.pop_front();
		visited[current] = true;

		//We push all unvisited children onto the stack
		for (uint_t i = 0; i < size; i++) {
			//Forward edges
			if (am[current][i] == 1 && visited[i] == false) {
				stack.push_front(i);
			}
			//Backward edges
			if (am[i][current] == 1 && visited[i] == false) {
				stack.push_front(i);
			}

		}
	}

	//Test to see if every node was visited or not	
	int connected = 0;
	for (uint_t it = 0; it < size; it++) {
		if (visited[it] == false) {
			connected = (int) it;
		}
	}
	return connected;
}

//Generates an adjacency matrix representation of a DAG
void DagGen::gen_erdos_matrix(matrix_t &am, uint_t size, double p) {

	//Get size of DAG to generate
	uint_t n = size;
	
	//Declare an adjacency matrix and initialize it to zeroes
	am.resize(n);
	for (uint_t i = 0; i < n; i++) {
		am[i].resize(n);
	}

	for (uint_t i = 0; i < n; i++) {
		for (uint_t j = 0; j < n; j++) {
			am[i][j] = 0;
		}
	}

	uint_t node;
	uint_t target;
	int graph_is_connected;

	for (uint_t i = 0; i < n; i++) {
		for (uint_t j = i + 1; j < n; j++) {
			//			double tmp = rand()/(RAND_MAX+1.0);
			// Son -- begin
			double tmp = Common::uniform_real_gen(0, 1);
			// Son -- end

			if(tmp <= p)
			{
				//Check to see if adding this edge adds a cycle
				//if(!acyclic_check(am,n,i,j)){
				//	cout << "cycle occurs" << endl; abort();}
				am[i][j] = 1;
			}
		}
	}
	graph_is_connected = weak_conn_add(am);	
	//	if(graph_is_connected != 0){cout<<".";}
	//We continue to add edges (i.e. loop) until we have a DAG that weakly
	//connects all vertices
	while (graph_is_connected != 0) {
		target = graph_is_connected;
		//		node = rand()%target;
		// Son -- begin
		node = Common::uniform_int_gen(0, target-1);
		// Son -- end

		//Check to see if adding this edge adds a cycle
		//if(acyclic_check(am,n,node,target)){
			//If no cycle is introduced, add the edge to the graph
			am[node][target] = 1;
			//If the graph is now weakly connected then we terminate
			graph_is_connected = weak_conn_add(am);
		//} else {
			//If a cycle is introduced, mark it as such
		//	cout << "cycle occurs" << endl; abort();
		//}	
	}

	/** For debugging **
		cout << "Resultant graph:" << endl;
		for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
		cout << am[i][j] << " ";
		}
		cout << endl;
		}
	*/
	
}

//Generates a set of random computational work values to go with 
//generated adjacency matrix to form a complete task DAG
//Return the total computational load of all nodes
uint_t DagGen::gen_node_lengths(vector<uint_t>& pa, uint_t size, uint_t min, uint_t max) {
	uint_t temp;
	uint_t C = 0;
	
	if (min == 0) {
		cout << "Minimum period length must be positive" << endl;
		abort();
	}

	pa.clear();
	pa.resize(size);

	for (uint_t i = 0; i < size; i++) {
		//		temp = min + rand()%(max-min+1);
		// Son -- begin
		temp = Common::uniform_int_gen(min, max);
		// Son -- end

		pa[i] = temp;
		C += temp;
	}
	
	return C;
}

uint_t DagGen::gen_periods_NP(vector<uint_t>& pa, uint_t size, uint_t min, uint_t max) {
	uint_t temp;
	uint_t C = 0;
	
	if (min == 0) {
		cout << "Minimum period length must be positive" << endl;
		abort();
	}

	pa.resize(size);

	int NP = (int)(max/min);
	for (uint_t i = 0; i < size; i++) {
		temp = min *(1 + rand()%NP);
		pa[i] = temp;
		C += temp;
	}

	return C;
}

//Calculates the span of a DAG
uint_t DagGen::calc_span_old(const matrix_t& am, const vector<uint_t>& pa) {

	uint_t size = pa.size();
	
	list<uint_t> stack;
	uint_t current;
	list<uint_t> source_nodes;
	vector<uint_t> dist (size, 0); 
	uint_t span = 0;	
	
	//Find all source nodes in the graph
	//If a column in the adjacency matrix is zeroes, then that node is a source
	for (uint_t j = 0; j < size; j++) {
		bool source = true;
		for (uint_t i = 0; i < size; i++) {
			if (am[i][j] == 1) {
				source = false;
			}
		}
		if (source) {
			source_nodes.push_back(j);
		}
	}
	
	//Initialize the starting nodes to be source nodes
	for (list<uint_t>::iterator iter = source_nodes.begin(); iter != source_nodes.end(); iter++) {
		int a_source = *iter;
		stack.push_back(a_source);
		dist[a_source] = pa[a_source];
	}
	
	uint_t temp;
	
	while (!stack.empty()) {
		current = stack.front();
		stack.pop_front();
		
		for (uint_t i = 0; i < size; i++) {
			if (am[current][i] == 1 ) {
				
				temp = dist[current] + pa[i];      
				
				//Push back the starting times of children nodes
				//This is not optimal because we process multiple
				//nodes more than once
				if (temp > dist[i]) {
					dist[i] = temp;
					stack.push_back(i);
					if (span < temp) {
						span = temp;
					}
				}
			}
		}
	}
	
	return span;
}

// Each pair in the output distance vector contains the longest distance to a node and its parent.
// For source nodes we use -1 for their parent nodes.
// INPUT: @am: the adjacency matrix of the DAG
// @pa: the WCETs of the subtasks.
// OUTPUT: @dist: vector of longest distances to subtasks from a source.
// @crit_path: a critical-path of the DAG.
// Return: the critical-path length of this DAG.
uint_t DagGen::calc_span(const matrix_t& am, const vector<uint_t>& pa,
						 vector<pair<uint_t, int>>& dist, vector<uint_t>& crit_path) {
	uint_t size = pa.size();

	dist.clear();
	dist.resize(size);
	// Initialize the distances and parents.
	for (uint_t i = 0; i < size; ++i) {
		dist[i].first = pa[i]; // longest distance.
		dist[i].second = -1; // parent.
	}

	// Initialize the span to be the longest wcet of individual subtasks.
	uint_t span = 0;
	for (uint_t wcet : pa) {
		span = max(wcet, span);
	}
	
	// Go through nodes in topological order.
	for (uint_t i = 0; i < size; ++i) {
		for (uint_t j = i + 1; j < size; ++j) {
			if (am[i][j] == 1) {
				// Node j is an adjacent of node i.
				if (dist[j].first < dist[i].first + pa[j]) {
					dist[j].first = dist[i].first + pa[j];
					dist[j].second = i;
					if (span < dist[j].first) {
						span = dist[j].first;
					}
				}
			}
		}
	}

	// Compute an arbitrary critical path of the DAG.
	stack<uint_t> temp;
	for (int i = 0; i < dist.size(); ++i) {
		if (dist[i].first == span) {
			uint_t curr = i;
			while (curr != -1) {
				temp.push(curr);
				curr = dist[curr].second;
			}
			break;
		}
	}

	// Convert to a vector for containing the critical path.
	crit_path.clear();
	while (!temp.empty()) {
		crit_path.push_back(temp.top());
		temp.pop();
	}
	
	return span;
}


// Return the total work of all nodes that are reachable from a given node.
uint_t DagGen::calc_subdag_work(uint_t node, const matrix_t& am, const vector<uint_t>& wcets) {
	uint_t sum_work = 0;
	unordered_set<uint_t> visited; // set of visited nodes.
	queue<uint_t> frontier;
	frontier.push(node);
	visited.insert(node);
	uint_t size = wcets.size();
	
	// Use BFS to find all reachable nodes.
	while (!frontier.empty()) {
		uint_t a = frontier.front();
		frontier.pop();
		sum_work += wcets[a];
		for (uint_t i = a + 1; i < size; ++i) {
			if (am[a][i] == 1 && visited.find(i) == visited.end()) {
				frontier.push(i);
				visited.insert(i);
			}
		}
	}

	return sum_work;
}

// Return a longest path from a given subtask to sink, together with its length.
pair<vector<uint_t>, uint_t> DagGen::calc_longest_path(uint_t sub_id, const matrix_t &am, const vector<uint_t> &wcets) {
	uint_t size = wcets.size();
	// For each subtask, store the longest path lengths from the given subtask. 
	// Its immediate predecessor is stored in the second value of the pair.
	vector<pair<uint_t, int>> path_lengths(size);
	for (int i = 0; i < size; ++i) {
		// Length of longest path from the given node.
		if (i == sub_id) {
			path_lengths[i].first = wcets[i];
		} else {
			path_lengths[i].first = 0;
		}
		path_lengths[i].second = -1; // Its parent node in the path. -1 means no parent.
	}
	// Length of longest path among all paths from the given node.
	uint_t max_length = wcets[sub_id];

	for (uint_t i = sub_id; i < size; ++i) {
		if (path_lengths[i].first == 0) {
			// Ignore subtasks that are not reachable from sub_id.
			continue;
		}
		
		for (uint_t j = i + 1; j < size; ++j) {
			if (am[i][j] == 1) { // There is an edge from i to j.
				if (path_lengths[j].first < path_lengths[i].first + wcets[j]) {
					path_lengths[j].first = path_lengths[i].first + wcets[j];
					path_lengths[j].second = i; // Set j's parent to be i.
					if (max_length < path_lengths[j].first) {
						max_length = path_lengths[j].first;
					}
				}
			}
		}
	}

	// Construct the longest path.
	stack<uint_t> temp;
	for (uint_t i = sub_id; i < size; ++i) {
		if (path_lengths[i].first == max_length) {
			uint_t curr = i;
			while (curr != -1) {
				temp.push(curr);
				curr = path_lengths[curr].second;
			}
			break;
		}
	}

	// Convert the path to vector and return.
	vector<uint_t> path;
	while (!temp.empty()) {
		path.push_back(temp.top());
		temp.pop();
	}

	return {path, max_length};
}

// For debugging the generated DAG.
void DagGen::test() {
	matrix_t am{{0, 0, 1, 1, 1, 1, 1, 1, 1, 1},
				{0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
				{0, 0, 0, 1, 1, 1, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
	vector<uint_t> wcets{3, 18, 18, 1, 1, 1, 1, 1, 16, 16};
	uint_t size = wcets.size();
	uint_t work = accumulate(wcets.begin(), wcets.end(), 0);
	vector<pair<uint_t, int>> dist;
	vector<uint_t> crit_path;
	uint_t span = calc_span(am, wcets, dist, crit_path);

	cout << "DAG Adjacency Matrix:" << endl;
	for (auto& vec : am) {
		for (uint_t a : vec) {
			cout << a << " ";
		}
		cout << endl;
	}

	cout << "Critical path: ";
	for (uint_t i : crit_path) {
		cout << "(" << i << ", " << wcets[i] << ") ";
	}
	cout << endl;
	
	cout << "DAG work: " << work << ", critical path length: " << span << endl;
	for (int i = 0; i < size; ++i) {
		cout << "Subdag work at node " << i << ": " << calc_subdag_work(i, am, wcets) << endl;
	}
}


//Generates the DAG used as an example in Preemptive and Non-Preemptive 
//Scheduling for Parallel Tasks of DAG Model by Abu and Kunal
//Note that nodes here are numbered one less than the paper
void DagGen::test_DAG(matrix_t& am, vector<uint_t>& pa){

	am.resize(10);
	for (uint_t i = 0; i < 10; i++) {
		am[i].resize(10);
		for (uint_t j = 0; j < 10; j++) {
			am[i][j] = 0;
		}
	}	
	
	pa.resize(10);
	
	//Specify outgoing edges
	am[0][3] = 1;
	am[0][6] = 1;
	am[0][7] = 1;
	am[0][9] = 1;
	
	am[1][3] = 1;
	am[1][6] = 1;
	am[1][4] = 1;
	
	am[2][5] = 1;
	am[2][8] = 1;

	am[3][7] = 1;

	am[4][6] = 1;
	am[4][5] = 1;

	am[5][8] = 1;

	am[6][7] = 1;

	am[7][8] = 1;
	am[7][9] = 1;

	//Specify periods
	pa[0] = 4;	
	pa[1] = 2;
	pa[2] = 4;
	pa[3] = 5;
	pa[4] = 3;
	pa[5] = 4;
	pa[6] = 2;
	pa[7] = 2;
	pa[8] = 3;
	pa[9] = 3;
}
