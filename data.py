# A python script to grab the data
# import sys
import os
import numpy as np
import pandas as pd
from scipy.io import loadmat


class NOData:

    def __init__(self):
        self._path = os.path.join(os.getcwd(), 'data')
        self._sessions, self._session_nrs = self._define_session()
        self._markers = self._define_event_markers()

    @property
    def path(self):
        return self._path

    @property
    def sessions(self):
        return self._sessions

    @property
    def session_nrs(self):
        return self._session_nrs

    @property
    def markers(self):
        return self._markers

    def _hack_mat_data_structure(self, session_line):
        """
        The original session data is saved as .mat file using the matlab program.
        read in .mat files using the scipy.io loadmat package has a special data structure.

        we want to use this program to hack in the data structure and turn it into a dictionary.
        """

        session = session_line[0][0][0]
        session_id = session_line[0][1][0]
        experiment_id_learn = session_line[0][2][0][0]
        experiment_id_recog = session_line[0][3][0][0]
        task_descr = session_line[0][4][0]
        variant = session_line[0][5][0][0]
        block_id_learn = session_line[0][6][0][0]
        block_id_recog = session_line[0][7][0][0]
        patient_nb = session_line[0][8][0][0]
        patient_session = session_line[0][9][0][0]
        diagnosis_code = session_line[0][10][0][0]
        return self._make_sess_dict(session, session_id, experiment_id_learn, experiment_id_recog,
                                    task_descr, variant, block_id_learn, block_id_recog, patient_nb,
                                    patient_session, diagnosis_code)

    def _construct_data_path(self, session_name, target_folder):
        """
            A method used to construct paths to desired data.
        """
        path = os.path.join(self.path, target_folder, session_name, 'NO')
        return path

    def _extract_event_periods(self, events, experiment_type='recog'):
        """
        A method to separate the important events from the raw event files.
        :param events: The raw event data, output from _get_event_data
        :param experiment_type: recog or learn experiment_type
        :return: a list of four dataframes for important events
        """
        # if experiment_type == 'recog':
        stimulus_on_events = events.loc[events[1] == self._markers['stimulus_on']]
        stimulus_off_events = events.loc[events[1] == self._markers['stimulus_off']]
        question_on_events = events.loc[events[1] == self._markers['delay1_off']]
        trial_end_events = events.loc[events[1] == self._markers['delay2_off']]
        recog_response_events = events.loc[(events[1] == self._markers['response_1']) |
                                           (events[1] == self._markers['response_2']) |
                                           (events[1] == self._markers['response_3']) |
                                           (events[1] == self._markers['response_4']) |
                                           (events[1] == self._markers['response_5']) |
                                           (events[1] == self._markers['response_6'])]

        if (stimulus_on_events.shape[0] != 100) | (stimulus_off_events.shape[0] != 100) | \
           (question_on_events.shape[0] != 100) | (trial_end_events.shape[0] != 100):
            print('Somethings wrong with the event data, each event should have 100 trials.')
            raise ValueError

        events_time_points = pd.DataFrame({'stimulus_on': np.asarray(stimulus_on_events[0]),
                                           'stimulus_off': np.asarray(stimulus_off_events[0]),
                                           'question_on': np.asarray(question_on_events[0]),
                                           'trial_end': np.asarray(trial_end_events[0]),
                                           # raw response from 31 - 36, designed response 1 - 6
                                           'recog_responses': np.asarray(recog_response_events[1] - 30)})
        # TODO experiment_type = 'learn'
        #print(events_time_points)
        return events_time_points

    def _get_event_data(self, session_nr, experiment_type='recog'):
        """
        Load event data from the raw data set with desired experiment ID
        """
        session_name = self.sessions[session_nr]['session']
        experiment_id_recog = self.sessions[session_nr]['experiment_id_recog']
        experiment_id_learn = self.sessions[session_nr]['experiment_id_learn']
        event_path = os.path.join(self._construct_data_path(session_name, 'events'), 'eventsRaw.mat')

        events = pd.DataFrame(loadmat(event_path)['events'])
        if experiment_type == 'recog':
            events = events.loc[events[2] == experiment_id_recog]
        elif experiment_type == 'learn':
            events = events.loc[events[2] == experiment_id_learn]
        else:
            print('please enter a correct option for experiment ID, now return the entire matrix')
        print(events.shape)
        print(events.tail())
        return events

    def get_trials_from_cell(self, cell, experiment_type='recog'):
        """
        This method takes in the cell raw spike data and separate it into different trials based on the event
        timestamps.
        :param cell: A cell object contains the cell information including the raw spike train.
        :param experiment_type: TODO: learning = 'learning' or recognition phase = 'recog'
        :return: A 100 x 7 list with columns indexes:
                 0 = category_recog_number
                 1 = category_recog_desc
                 2 = img_file_path
                 3 = stimuli_id
                 4 = new_old_recog
                 5 = Response dataframe for scale from 1 - 6 for the recognition phase,
                 6 = a list contains the spikes for each trials, 100 trials.
        """

        # get labels from cell
        trial_label_recog = self._get_labels_from_cell(cell)

        session_nr = cell.session_nr
        events = self._get_event_data(session_nr, experiment_type=experiment_type)
        events_time_points = self._extract_event_periods(events, experiment_type=experiment_type)
        recog_responses = events_time_points['recog_responses']

        cell_raw_spike_timestamps = pd.DataFrame(cell.raw_spike_timestamps)
        for i in range(0, 100):
            trial_start = events_time_points.iloc[i, 0]
            trial_end = events_time_points.iloc[i, 3]
            trial = cell_raw_spike_timestamps.loc[(cell_raw_spike_timestamps[0] > trial_start) &
                                                  (cell_raw_spike_timestamps[0] < trial_end)]
            trial = trial[0].tolist()
            trial_label_recog[i].append(recog_responses[i])
            trial_label_recog[i].append(np.asarray(trial) - trial_start)
        return trial_label_recog

    def _get_labels_from_cell(self, cell, experiment_type='recog'):
        """
        This method use the original experiment stimuli data set to extract the order and categories of the stimuli
        shown to each subjects. The label of the stimulus for the new old recognition task is also obtained using
        this method.
        :param cell: a cell object contains the session number information to extract necessary data from the session
                     dictionary.
        :param experiment_type: 'recog' = recognition phase, 'learn' = learning phase
        :return: a 100 x 5 list with columns 0 = category_recog_number, 1 = category_recog_desc, 2 = img_file_path,
                 3 = stimuli_id, 4 = new_old_recog
        """
        session_nr = cell.session_nr
        block_id_recog = self.sessions[session_nr]['block_id_recog']
        block_id_learn = self.sessions[session_nr]['block_id_learn']
        variant = self.sessions[session_nr]['variant']
        filename = 's'
        if variant == 1:
            filename = 'NewOldDelay_v3.mat'
            filename2 = 'NewOldDelayStimuli.mat'
        elif variant == 2:
            filename = 'NewOldDelay2_v3.mat'
            filename2 = 'NewOldDelayStimuli2.mat'
        elif variant == 3:
            filename = 'NewOldDelay3_v3.mat'
            filename2 = 'NewOldDelayStimuli3.mat'

        path_to_labels = os.path.join(self.path, 'stimFiles', filename)
        experiment_stimuli = loadmat(path_to_labels)['experimentStimuli']
        stimuli_recog = experiment_stimuli[0, block_id_recog-1][3]
        new_old_recog = experiment_stimuli[0, block_id_recog-1][4]
        stimuli_learn = experiment_stimuli[0, block_id_learn-1][2]

        path_to_categories = os.path.join(self.path, 'stimFiles', filename2)
        category_mat = loadmat(path_to_categories)
        category = category_mat['categories']
        category_mapping = pd.DataFrame(category_mat['categoryMapping']).set_index(0)
        file_mapping = category_mat['fileMapping']

        category_file_label_recog = []
        for i in range(0, 100):
            temp = []
            stimuli_recog_value = stimuli_recog[0][i]
            category_number = category_mapping.loc[stimuli_recog_value, 1]
            temp.append(category_number)
            temp.append(category[0, category_number-1][0])
            file_path = file_mapping[0][i][0].replace('C:\code\\', '')
            temp.append(file_path)
            temp.append(stimuli_recog_value)
            temp.append(new_old_recog[0][i])

            category_file_label_recog.append(temp)



        # print(category.shape)
        # print(category_mapping)
        # print(file_mapping.shape)

        return category_file_label_recog

    def _define_session(self):
        """ 
        Create a dictionary to contain the usable session information.
        """

        # read in session info
        session_path = os.path.join(self.path, 'sessions.mat')
        sessions_mat = loadmat(session_path)
        n = sessions_mat['NOsessions'].shape[1]

        # declare variables to record session data
        sessions = {}
        session_nrs = []

        for i in range(0, n):
            session_line = sessions_mat['NOsessions'][:, i]
            if np.size(session_line[0][0]) != 0:
                sessions[i] = self._hack_mat_data_structure(session_line)
                session_nrs.append(i)
        return sessions, session_nrs

    def ls_cells(self, session_nr):
        """
        The ls_cells function list all the available cells for a particular session number
        input:
            session_nr = the session number in the dictionary keys
        output:
            cell_list = a list of tuples (channel_nr, cluster_id)
        """

        session_name = self.sessions[session_nr]['session']
        brain_area_file_path = os.path.join(self._construct_data_path(session_name, 'events'), 'brainArea.mat')
        brain_area = loadmat(brain_area_file_path)['brainArea']
        cell_list = []

        for i in range(0, brain_area.shape[0]):
            cell_list.append((brain_area[i][0], brain_area[i][1]))
        return cell_list

    def pop_cell(self, session_nr, channelnr_clusterid):
        """
        This method pops a particular cell to the user
        input:
            session_nr = the session number that we would like to use
            channelnr_clusterid = the tuple contains both the channel id and the cluster id to select the
                                    desired cell
        output:
            cell = a cell object
        """

        session_name = self.sessions[session_nr]['session']
        brain_area_file_path = os.path.join(self._construct_data_path(session_name, 'events'), 'brainArea.mat')
        brain_area = loadmat(brain_area_file_path)['brainArea']
        df_brain_area = pd.DataFrame(brain_area)

        brain_area_cell = df_brain_area.loc[
            (df_brain_area[0] == channelnr_clusterid[0]) & (df_brain_area[1] == channelnr_clusterid[1])]

        cell_path = os.path.join(self._construct_data_path(session_name, 'sorted'),
                                 'A' + str(channelnr_clusterid[0]) + '_cells.mat')
        cell = Cell(cell_path, session_nr, session_name, *np.asarray(brain_area_cell))
        return cell

    def test(self):
        cell = self.pop_cell(114, (17, 1))

        #trials, recog_response = self.get_trials_from_cell(cell)

        stimuli_recog, new_old_recog = self.get_labels_from_cell(cell)

        # print(recog_response)
        # print(stimuli_recog)
        # print(new_old_recog)

    # Static helper methods
    @staticmethod
    def _make_sess_dict(session, session_id, experiment_id_learn, experiment_id_recog,
                        task_descr, variant, block_id_learn, block_id_recog, patient_nb,
                        patient_session, diagnosis_code):
        """
        :param session: The session name
        :param session_id: The session id
        :param experiment_id_learn: The experiment id that is for the learning task
        :param experiment_id_recog: The experiment id that is for recognition task
        :param task_descr:
        :param variant: The set of images used for the experiment
        :param block_id_learn: The block id used to index the labels for learning
        :param block_id_recog: The block id used to index the labels for recognition
        :param patient_nb: The patient id
        :param patient_session: The session where the patient is at
        :param diagnosis_code: The code that shows whether a patient has epilepsy or not
        :return: A session dictionary that contains all the above session information
        """
        session = {'session': session,
                   'session_id': session_id,
                   'experiment_id_learn': experiment_id_learn,
                   'experiment_id_recog': experiment_id_recog,
                   'task_descr': task_descr,
                   'variant': variant,
                   'block_id_learn': block_id_learn,
                   'block_id_recog': block_id_recog,
                   'patient_nb': patient_nb,
                   'patient_session': patient_session,
                   'diagnosis_code': diagnosis_code}

        return session

    @staticmethod
    def _define_event_markers():
        """
        This static method is used to define the useful markers to index the event data.
        :return: marker dictionary
        """
        markers = {'stimulus_on': 1,
                   'stimulus_off': 2,
                   'delay1_off': 3,
                   'delay2_off': 6,
                   'response_1': 31,
                   'response_2': 32,
                   'response_3': 33,
                   'response_4': 34,
                   'response_5': 35,
                   'response_6': 36,
                   'response_learning_animal': 21,
                   'response_learning_non_animal': 22,
                   'experiment_on': 55,
                   'experiment_off': 66}
        return markers


class Cell:

    def __init__(self, cell_path, session_nr, session_name, cell_info):
        self._data_path = cell_path
        self._session_nr = session_nr
        self._session_name = session_name
        self._channel_nr = cell_info[0]
        self._cluster_id = cell_info[1]
        self._ori_cluster_id = cell_info[2]
        self._brain_area_id = cell_info[3]
        self._raw_spike_timestamps = self._load_cell_data(cell_path)

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

    def _load_cell_data(self, cell_path):
        """
        load the raw cell data and capture the spike train timestamps
        """
        cell_data = loadmat(cell_path)
        channel_raw_spike_timestamps = pd.DataFrame(cell_data['spikes'])
        cell_raw_spike_timestamps = np.asarray(
            channel_raw_spike_timestamps.loc[channel_raw_spike_timestamps[0] == self.cluster_id, 2])
        return cell_raw_spike_timestamps






