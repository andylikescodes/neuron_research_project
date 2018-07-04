# The cell.py file contains the cell associated classes
from scipy.stats import f_oneway
import numpy as np

class Cell:

    def __init__(self, cell_path, session_nr, session_name, cell_info, raw_spike_timestamps, trials):
        self._data_path = cell_path
        self._session_nr = session_nr
        self._session_name = session_name
        self._channel_nr = cell_info[0]
        self._cluster_id = cell_info[1]
        self._ori_cluster_id = cell_info[2]
        self._brain_area_id = cell_info[3]
        self._raw_spike_timestamps = raw_spike_timestamps
        self._trials = trials

    @property
    def data_path(self):
        return self._data_path

    @property
    def session_nr(self):
        return self._session_nr

    @property
    def session_name(self):
        return self._session_name

    @property
    def channel_nr(self):
        return self._channel_nr

    @property
    def cluster_id(self):
        return self._cluster_id

    @property
    def ori_cluster_id(self):
        return self._ori_cluster_id

    @property
    def brain_area_id(self):
        return self._brain_area_id

    @property
    def raw_spike_timestamps(self):
        return self._raw_spike_timestamps

    @property
    def trials(self):
        return self._trials

    def vs_test(self):
        """
        Anova test for stimuli categories
        :return: p value for the anova test
        """
        trials = (trial for trial in self.trials)
        cat_1 = []
        cat_2 = []
        cat_3 = []
        cat_4 = []
        cat_5 = []
        for trial in trials:
            if trial.category == 1:
                cat_1.append(trial.win_spike_rate(200, 1700))
            elif trial.category == 2:
                cat_2.append(trial.win_spike_rate(200, 1700))
            elif trial.category == 3:
                cat_3.append(trial.win_spike_rate(200, 1700))
            elif trial.category == 4:
                cat_4.append(trial.win_spike_rate(200, 1700))
            elif trial.category == 5:
                cat_5.append(trial.win_spike_rate(200, 1700))
        return f_oneway(cat_1, cat_2, cat_3, cat_4, cat_5)[1]

    def baseline_test(self):
        """
        Anova test to test the difference between the baseline rate and the stim on rate
        :return: p value of the anova test
        """
        trials = (trial for trial in self.trials)
        baseline = []
        stim_period = []
        for trial in trials:
            baseline.append(trial.baseline_spike_rate())
            stim_period.append(trial.win_spike_rate(200, 1700))
        return f_oneway(baseline, stim_period)[1]

    def ms_test(self, n):
        """
        A bootstrap test for new old test
        :param n: number of bootstraps
        :return: a p value of the bootstrap test
        """
        trials = (trial for trial in self.trials)

        old = np.array([])
        new = np.array([])

        for trial in trials:
            if (trial.new_old_recog == 0) & (trial.response_recog <= 3.):
                new = np.append(new, trial.win_spike_rate(200, 1700))
            elif (trial.new_old_recog == 1) & (trial.response_recog >= 4.):
                old = np.append(old, trial.win_spike_rate(200, 1700))

        m_ = len(new)
        n_ = len(old)

        all_m = np.mean(np.concatenate([new, old]))

        new_m = new - np.mean(new) + all_m
        old_m = old - np.mean(old) + all_m
        new_bootstrap = np.array([])
        old_bootstrap = np.array([])
        for i in range(0, n):
            random_ints = np.random.randint(m_, size=m_)
            new_samples = new_m[random_ints]
            new_bootstrap = np.append(new_bootstrap, np.mean(new_samples))

            random_ints = np.random.randint(n_, size=n_)
            old_samples = old_m[random_ints]
            old_bootstrap = np.append(old_bootstrap, np.mean(old_samples))

        t = np.abs(new_bootstrap - old_bootstrap)
        t_obs = np.abs(np.mean(new) - np.mean(old))

        return np.mean(t >= t_obs)

