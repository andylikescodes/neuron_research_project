# Trial class


class Trial:
    def __init__(self, category, new_old_recog, response_recog, category_name,
                 file_path, stimuli_id, trial_duration, trial_timestamps):
        self._category = category
        self._new_old_recog = new_old_recog
        self._response_recog = response_recog
        self._category_name = category_name
        self._file_path = file_path
        self._stimuli_id = stimuli_id
        self._trial_duration = trial_duration
        self._trial_timestamps = trial_timestamps
        self._total_spike_counts = len(self._trial_timestamps)

    @property
    def category(self):
        return self._category

    @property
    def new_old_recog(self):
        return self._new_old_recog

    @property
    def response_recog(self):
        return self._response_recog

    @property
    def category_name(self):
        return self._category_name

    @property
    def file_path(self):
        return self._file_path

    @property
    def stimuli_id(self):
        return self._stimuli_id

    @property
    def trial_duration(self):
        return self._trial_duration

    @property
    def trial_timestamps(self):
        return self._trial_timestamps

    def win_spike_count(self, win_start, win_end):
        """
        Calculate the spike count in a window
        :param win_start: window starting time in millisecond
        :param win_end: window ending time in millisecond
        :return: spike count in the window
        """
        start = win_start * 1000
        end = win_end * 1000

        timestamps_within_window = self._trial_timestamps[(self._trial_timestamps > start) *
                                                          (self._trial_timestamps < end)]
        return len(timestamps_within_window)

    def win_spike_rate(self, win_start, win_end):
        """
        Calculate the spike count rate in a window
        :param win_start:
        :param win_end:
        :return: spike count rate in the window
        """
        start = win_start * 1000
        end = win_end * 1000

        timestamps_within_window = self._trial_timestamps[(self._trial_timestamps > start) *
                                                          (self._trial_timestamps < end)]
        return len(timestamps_within_window) / ((end - start) / 1000000)
