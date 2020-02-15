HOST=$(shell hostname)
ifeq ($(HOST),shakespeare.cse.wustl.edu)
GXX=/usr/local/gcc5/bin/g++
else
GXX=g++
endif

sched: *.cpp
	$(GXX) -g -std=c++11 main.cpp dag_task.cpp rand_dag.cpp proposed_scheduler.cpp naive_greedy_scheduler.cpp \
	baruah_greedy_scheduler.cpp taskset.cpp federated_common.cpp semi_federated.cpp \
	reservation_federated.cpp baruah_federated.cpp proposed_federated.cpp stretching_algorithm.cpp \
	decomp_algorithm.cpp common.cpp -o sched

clean:
	rm sched
