# A script for testing the functions
# Now it includes a function to plot the raster plots of cells
# Cells can be grab from the data set anytime. Just enter the session
# id and the channel and cluster id for the cells

# It will be continuely refined and a more organized class for analysis
# will be added later.

from data import NOData
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

NO_data = NOData()
cell = NO_data.pop_cell(114, (17, 1))
trial_labels = NO_data.get_trials_from_cell(cell)

print(trial_labels[0])

colors1 = np.array([[1, 0, 0],
                    [0, 1, 1],
                    [0, 0, 1],
                    [0, 1, 0],
                    [1, 0, 1]])

sorted_new_data = sorted(trial_labels, key=itemgetter(0))


def grab_trials(data, colors):
    trials_data = []
    colors_mapping = []
    for i in range(0, 100):
        trials_data.append(data[i][6])
        colors_mapping.append(colors[data[i][0]-1])
    return trials_data, colors_mapping


trials, colors = grab_trials(sorted_new_data, colors1)

print(trials[0])
print(colors[0])

plt.eventplot(np.asarray(trials)/1000, colors=colors)

plt.axvspan(1000, 2000, color='grey', alpha=0.2)

plt.show()

