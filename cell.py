# The cell.py file contains the cell associated classes


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


