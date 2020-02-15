from tkinter import *
import sys


# Show the adjacency list for the DAG task together with its parameters.
class AdjListDrawing(Frame):
    def __init__(self, num_subtasks, adj_list, wcets, span, deadline, period, master=None):
        super().__init__(master)
        self.num_subtasks = num_subtasks
        self.adj_list = adj_list
        self.wcets = wcets
        self.span = span
        self.deadline = deadline
        self.period = period
        self.pack()
        self.show()

    def show(self):
        # 3 heading lines plus 1 line for each subtask.
        num_lines = 3 + int(self.num_subtasks)
        self.app = Text(self, height=num_lines) # Default is 24 lines.
        self.app.pack()

        self.app.insert(END, "Number of subtasks: " + str(self.num_subtasks) + "\nExecution times: ")
        self.work = 0
        wcets_line = ""
        for node_len in self.wcets:
            self.work += int(node_len)
            wcets_line += node_len + " "
        wcets_line += "\n"
        para_line = "Work: " + str(self.work) + ". Span: " + str(self.span) + ". Deadline: " + \
        str(self.deadline) + ". Period: " + str(self.period) + "\n"
        self.app.insert(END, wcets_line)
        self.app.insert(END, para_line)

        for sub_id, adj in self.adj_list.items():
            line = "v" + str(sub_id) + " ==> ["
            for neighbor in adj:
                line += "v" + neighbor + " "
            line += "]\n"
            self.app.insert(END, line)
        self.app.config(state=DISABLED)
        # Add label to this text widget.
        Label(self, text="DAG Structure", font=("Helvetica", 20)).pack()
        self.app.pack()

        
# Class for drawing the schedules rendered by different algorithms.        
class ScheduleDrawing(Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        # Row number for the schedules. This is used to call grid() method.
        self.row = 0

        
    # This and other drawing methods in this class are meant to be called by outside code.
    def draw_proposed_sched(self, num_subtasks, num_procs, subtask_to_frags, deadline):
        self.num_subtasks = num_subtasks
        self.num_procs = num_procs
        self.subtask_to_frags = subtask_to_frags
        self.deadline = deadline
        # Scale out the schedule.
        self.scale = 1
        if 1000 > int(self.deadline):
            self.scale = 1000//int(self.deadline)
        self.width = self.scale * int(self.deadline) + 100
        # The height of the proposed schedule.
        self.height = (self.num_procs + 1) * 50
        # Each schedule is put within a canvas.
        self.canvas = Canvas(self, width=self.width, height=self.height)
        self.canvas.grid(row=self.row, column=0)
        self.canvas.pack()
        self.row += 1 # Compute the row number for the next canvas.
        # Gap between lines.
        gap = 50
        x_start = 40
        x_end = self.width - 10
        # y-coordinate of the bottom line.
        bottom = self.height - 20
        lines = list()
        # Each line corresponds to a processor's timeline.
        for i in range(0, self.num_procs):
            y_coor = bottom - i * gap
            line_id = self.canvas.create_line(x_start, y_coor, x_end, y_coor, arrow=LAST)
            lines.append(line_id)
            label = "P" + str(i)
            self.canvas.create_text(x_start-20, y_coor-15, text=label)

        # Add a line for the DAG task's deadline.
        x_dln = x_start + self.scale * int(self.deadline)
        self.canvas.create_line(x_dln, 10, x_dln, bottom, arrow=LAST)
            
        # Draw schedule for the subtasks.
        rect_height = gap - 15
        for sub_id, fragments in self.subtask_to_frags.items():
            frag_no = 0
            for frag in fragments:
                proc_id = int(frag[0])
                scaled_start_time = int(frag[1]) * self.scale + x_start
                scaled_exe_time = int(frag[2]) * self.scale
                scaled_deadline = int(frag[3]) * self.scale
                # Top left coordinate.
                x_1 = scaled_start_time
                y_1 = bottom - proc_id * gap - rect_height
                # Bottom right coordinate.
                x_2 = x_1 + scaled_exe_time
                y_2 = y_1 + rect_height
                self.canvas.create_rectangle(x_1, y_1, x_2, y_2, fill='cyan')
                # Add label to the rectangles.
                label = "v" + str(sub_id) + "_" + str(frag_no)
                frag_no += 1
                self.canvas.create_text((x_1+x_2)/2, (y_1+y_2)/2, text=label, width=scaled_exe_time)
                # Add time ticks to the bottom line.
                if proc_id == 0:
                    self.canvas.create_text(scaled_start_time, bottom + 15, text=str(frag[1]))
                    if int(frag[1]) + int(frag[2]) == int(self.deadline):
                        self.canvas.create_text(scaled_start_time + scaled_exe_time, bottom + 15, \
                                           text=self.deadline)
        # Add label to this frame.
        Label(self, text="Proposed Schedule", font=("Helvetica", 20)).pack()


    # Common method to draw both the naive schedule and the improved schedule by Baruah.
    # NOTE: call draw_proposed_sched() before this method since it needs some instance
    # variables to be initialized.
    def draw_greedy_sched(self, proc_to_subtasks, algo="baruah"):
        self.num_procs_greedy = len(proc_to_subtasks)
        greedy_height = (self.num_procs_greedy + 1) * 50
        canvas = Canvas(self, width=self.width, height=greedy_height)
        canvas.grid(row=self.row, column=0)
        canvas.pack()
        self.row += 1
        bottom = greedy_height - 20
        gap = 50
        x_start = 40
        x_end = self.width - 10
        rect_height = gap - 15
        # Add a line for the DAG task's deadline.
        x_dln = x_start + self.scale * int(self.deadline)
        canvas.create_line(x_dln, 10, x_dln, bottom, arrow=LAST)
        canvas.create_text(x_dln, bottom + 15, text=self.deadline)
        # Set of time ticks.
        tick_set = set()
        # Draw the rectangles for the subtask executions.
        for i in range(0, self.num_procs_greedy):
            y_coor = bottom - i * gap
            line_id = canvas.create_line(x_start, y_coor, x_end, y_coor, arrow=LAST)
            label = "P" + str(i)
            canvas.create_text(x_start-20, y_coor-15, text=label)
            subtasks = proc_to_subtasks[i]
            for subtask in subtasks:
                sub_id = int(subtask[0])
                scaled_start_time = self.scale * int(subtask[1]) + x_start
                scaled_wcet = self.scale * int(subtask[2])
                x_1 = scaled_start_time
                y_1 = bottom - i * gap - rect_height
                x_2 = x_1 + scaled_wcet
                y_2 = y_1 + rect_height
                canvas.create_rectangle(x_1, y_1, x_2, y_2, fill="cyan")
                # Add label to the rectangles.
                label = "v" + str(sub_id)
                canvas.create_text((x_1+x_2)/2, (y_1+y_2)/2, text=label, width=scaled_wcet)
                # Add time ticks to the bottom line (i.e., processor 0's line).
                if int(subtask[1]) not in tick_set:
                    canvas.create_text(scaled_start_time, bottom + 15, text=subtask[1])
                    tick_set.add(int(subtask[1]))
                end_time = int(subtask[1]) + int(subtask[2])
                if end_time not in tick_set:
                    canvas.create_text(scaled_start_time + scaled_wcet, bottom + 15, \
                                       text=str(end_time))
                    tick_set.add(end_time)

        if algo == "baruah":
            Label(self, text="Baruah Greedy Schedule", font=("Helvetica", 20)).pack()
        else:
            Label(self, text="Naive Greedy Schedule", font=("Helvetica", 20)).pack()


# The GUI's main class.
class Entry:
    # @file_path: path to the file containing the DAG structure.
    def __init__(self, dag_path, proposed_path, baruah_path=None, naive_path=None):
        # Read input files.
        self.read_dag(dag_path)
        self.read_proposed_sched(proposed_path)
        self.has_baruah_sched = False
        self.has_naive_sched = False
        # Map from processor to subtasks scheduled on this processor.
        self.proc_to_subtasks_baruah = dict()
        self.proc_to_subtasks_naive = dict()
        if baruah_path != None:
            self.has_baruah_sched = True
            self.read_greedy_sched(baruah_path, algo="baruah")
        if naive_path != None:
            self.has_naive_sched = True
            self.read_greedy_sched(naive_path, algo="naive")
            
        # Draw the DAG and schedules.
        self.run()

    def read_dag(self, dag_path):
        dag_f = open(dag_path, 'r')
        # First line: <number subtasks>
        self.num_subtasks = int(dag_f.readline())
        # Second line: <wcets>
        self.wcets = dag_f.readline().split()
        # Third line: <span> <deadline> <period>
        parameters = dag_f.readline().split()
        self.span = parameters[0]
        self.deadline = parameters[1]
        self.period = parameters[2]

        # Read adjacency list.
        self.adj_list = dict()
        for i in range(0, self.num_subtasks):
            line = dag_f.readline()
            tokens = line.split()
            sub_id = tokens[0]
            tokens.pop(0)
            self.adj_list[sub_id] = tokens
        dag_f.close()        

        
    def read_proposed_sched(self, proposed_path):
        f = open(proposed_path, 'r')
        # First line: <Number of subtasks>
        self.num_subtasks = int(f.readline())
        # Map from subtask id to a list of its fragments.
        self.subtask_to_frags = dict()
        self.num_procs = 1
        # Subsequent lines: schedule for each subtask.
        for i in range(0, self.num_subtasks):
            tokens = f.readline().split()
            num_frags = int(tokens[1])
            # List of fragments for this subtask, each is a tuple.
            frags = list()
            for j in range(0, num_frags):
                proc_id = tokens[2 + j * 4]
                self.num_procs = max(self.num_procs, int(proc_id) + 1)
                start_time = tokens[3 + j * 4]
                exe_time = tokens[4 + j * 4]
                deadline = tokens[5 + j * 4]
                frags.append((proc_id, start_time, exe_time, deadline))
            self.subtask_to_frags[i] = frags
        f.close()

    # Common method for reading greedy schedules.
    # Return a map from each processor to a list of subtasks scheduled on it.
    def read_greedy_sched(self, file_path, algo="baruah"):
        f = open(file_path, 'r')
        # First line: <Number of subtasks>
        no_procs = int(f.readline())
        # Subsequent lines: <processor id> <subtask id> <start time> <execution time>
        for i in range(0, no_procs):
            tokens = f.readline().split()
            no_subs = (len(tokens) - 1)//3;
            subtasks = list()
            for j in range(0, no_subs):
                sub_id = tokens[1 + 3 * j]
                start_time = tokens[2 + 3 * j]
                exe_time = tokens[3 + 3 * j]
                subtasks.append((sub_id, start_time, exe_time))
            if algo == "baruah":
                self.proc_to_subtasks_baruah[i] = subtasks
            else:
                self.proc_to_subtasks_naive[i] = subtasks
        
        
    def run(self):
        # A frame for showing the DAG structure.
        root = Tk()
        root.title("DAG Structure & Schedules")
        # A Frame for DAG's adjacency list.
        self.dag = AdjListDrawing(self.num_subtasks, self.adj_list, self.wcets, \
                                  self.span, self.deadline, self.period, master=root)

        # A Frame for showing the schedules.
        self.scheds = ScheduleDrawing(master=root)
        #outer_canvas = Canvas(root)
        #self.scheds = ScheduleDrawing(master=outer_canvas)
        self.scheds.draw_proposed_sched(self.num_subtasks, self.num_procs, self.subtask_to_frags, self.deadline)
        if self.has_baruah_sched == True:
            self.scheds.draw_greedy_sched(self.proc_to_subtasks_baruah, algo="baruah")
        if self.has_naive_sched == True:
            self.scheds.draw_greedy_sched(self.proc_to_subtasks_naive, algo="naive")

        # A scrollbar for viewing the schedules.
        #scrollbar = Scrollbar(root, orient="vertical", command=outer_canvas.yview)
        #outer_canvas.configure(yscrollcommand=scrollbar.set)
        #scrollbar.pack(side="right")
        #outer_canvas.pack()
        root.mainloop()
        

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print("ERROR: Usage: <program name> <DAG file> <proposed schedule file> [Baruah schedule file] [naive schedule file] !!\n")
        exit()
    dag_path = sys.argv[1]
    proposed_path = sys.argv[2]
    baruah_path = None
    naive_path = None
    if len(sys.argv) >= 4:
        baruah_path = sys.argv[3]
    if len(sys.argv) == 5:
        naive_path = sys.argv[4]
    entry = Entry(dag_path, proposed_path, baruah_path, naive_path)
