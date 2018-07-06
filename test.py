# A script for testing the functions
# Now it includes a function to plot the raster plots of cells
# Cells can be grab from the data set anytime. Just enter the session
# id and the channel and cluster id for the cells

# It will be continuously refined and a more organized class for analysis
# will be added later.

from data import NOData
import pickle
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import matplotlib
plt.switch_backend('module://backend_interagg')

NO_data = NOData()

# cell = NO_data.pop_cell(58, (9, 1))

f_all = open('data/all_cells.p', 'rb')
all_cells = pickle.load(f_all)

all_cells[0].psth('visual')
all_cells[0].raster()

# for session_nr in NO_data.session_nrs:
#     print('========> session number: ' + str(session_nr))
#     print('===> session name' + NO_data.sessions[session_nr]['session'])
#     cell_list = NO_data.ls_cells(session_nr)
#     for cell_id in NO_data.ls_cells(session_nr):
#         print('=> cell id: (' + str(cell_id[0]) + ', ' + str(cell_id[1]) + ')')
#         cell = NO_data.pop_cell(session_nr, cell_id)
#         if cell:
#             all_cell.append(cell)

# f = open('vs_cells.p', 'rb')
# vs_cells = pickle.load(f)
#
# f_ms = open('ms_cells.p', 'rb')
# ms_cells = pickle.load(f_ms)

# selected_cells = []
#
# vs_cells[0].psth()
# vs_cells[0].raster()

# for i in range(0, len(vs_cells)):
#     if vs_cells[i].vs_test() <= 0.000005:
#         selected_cells.append(vs_cells[i])
#
# for i in range(0, len(selected_cells)):
#     selected_cells[i].psth()
#     selected_cells[i].raster()


# vs_cells = []
#
# n_total = 0


# 98(4,1) 85(3,1) 58(9,1) somethings wrong with these
# trials = NO_data.get_trials_from_cell(cell)

# print(trials)
#
# colors1 = np.array([[1, 0, 0],
#                     [0, 1, 1],
#                     [0, 0, 1],
#                     [0, 1, 0],
#                     [1, 0, 1]])
#
# def grab_trials(data, colors):
#     trials_data = []
#     colors_mapping = []
#     for i in range(0, 100):
#         trials_data.append(data[i][6])
#         colors_mapping.append(colors[data[i][0]-1])
#     return trials_data, colors_mapping
#
#
# trials, colors = grab_trials(sorted_new_data, colors1)
#
# print(trials[0])
# print(colors[0])
#
# plt.eventplot(np.asarray(trials)/1000, colors=colors)
#
# plt.axvspan(1000, 2000, color='grey', alpha=0.2)
#
# plt.show()

