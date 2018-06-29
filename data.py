# A python script to grab the data
import sys, os
import numpy as np
import pandas as pd
from scipy.io import loadmat 


class NOData():

	def __init__(self):
		self.path = os.path.join(os.getcwd(), 'data')
		self.sessions, self.sessionNrs = self._define_session()

	# Private Methods:
	def _make_sess_dict(self, session, sessionID, EXPERIMENTIDLearn, EXPERIMENTIDRecog, \
						taskDescr, variant, blockIDLearn, blockIDRecog, patient_nb,\
						patientsession, diagnosisCode):
		session = {'session': session, \
					'sessionID': sessionID, \
					'EXPERIMENTIDLearn': EXPERIMENTIDLearn, \
					'EXPERIMENTIDRecog': EXPERIMENTIDRecog, \
					'taskDescr': taskDescr, \
					'variant': variant, \
					'blockIDLearn': blockIDLearn, \
					'patient_nb': patient_nb, \
					'patientsession': patientsession,\
					'diagnosisCode': diagnosisCode}

		return session

	def _hack_mat_data_structure(self, session_line):
		'''
			The original session data is saved as .mat file using the matlab program.
			read in .mat files using the scipy.io loadmat package has a special data structure.

			we want to use this program to hack in the data structure and turn it into a dictionary.
		'''

		session = session_line[0][0][0]
		sessionID = session_line[0][1][0]
		EXPERIMENTIDLearn = session_line[0][2][0][0]
		EXPERIMENTIDRecog = session_line[0][3][0][0]
		taskDescr = session_line[0][4][0]
		variant = session_line[0][5][0][0]
		blockIDLearn = session_line[0][6][0][0]
		blockIDRecog = session_line[0][7][0][0]
		patient_nb = session_line[0][8][0][0]
		patientsession = session_line[0][9][0][0]
		diagnosisCode = session_line[0][10][0][0]
		return self._make_sess_dict(session, sessionID, EXPERIMENTIDLearn, EXPERIMENTIDRecog, \
						taskDescr, variant, blockIDLearn, blockIDRecog, patient_nb,\
						patientsession, diagnosisCode)

	def _construct_data_path(self, session_name, target_folder):
		'''
			A method used to construct paths to desired data.
		'''
		path = os.path.join(self.path, target_folder, session_name, 'NO')
		return path

	def _get_event_data(self, sessionNr, experiment_type = 'recog'):
		'''
			Load event data from the raw dataset with desired experiment ID
		'''
		session_name = self.sessions[sessionNr]['session']
		EXPERIMENTIDRecog = self.sessions[sessionNr]['EXPERIMENTIDRecog']
		EXPERIMENTIDLearn = self.sessions[sessionNr]['EXPERIMENTIDLearn']
		event_path = os.path.join(self._construct_data_path(session_name, 'events'), 'eventsRaw.mat')
		
		events = pd.DataFrame(loadmat(event_path)['events'])
		print(events.shape)
		if experiment_type == 'recog':
			events = events.loc[events[2] == EXPERIMENTIDRecog]
		elif experiment_type == 'learn':
			events = events.loc[events[2] == EXPERIMENTIDLearn]
		else:
			print ('please enter a correct option for experiment ID, now return the entire matrix')
		print(events.shape)
		print(events.tail())
		return events

	def _define_session(self):
		'''
			Create a dictionary to contain the usable session information.
		'''
		
		# read in session info
		session_path = os.path.join(self.path, 'sessions.mat')
		sessions_mat = loadmat(session_path)
		n = sessions_mat['NOsessions'].shape[1]

		# declare variables to record session data
		sessions = {}
		sessionNrs = []

		for i in range(0,n):
			session_line = sessions_mat['NOsessions'][:,i]
			if np.size(session_line[0][0]) != 0:
				sessions[i] = self._hack_mat_data_structure(session_line)
				sessionNrs.append(i)
		return sessions, sessionNrs


	def ls_cells(self, sessionNr):
		'''
			The ls_cells fucntion list all the available cells for a particular session number
			input:
				sessionNr = the session number in the dictonary keys
			output:
				cell_list = a list of tuples (channelNr, clusterID)
		'''


		session_name = self.sessions[sessionNr]['session']
		print(session_name)
		brain_area_file_path = os.path.join(self._construct_data_path(session_name, 'events'), 'brainArea.mat')

		brainArea = loadmat(brain_area_file_path)['brainArea']

		cell_list = []

		for i in range(0, brainArea.shape[0]):
			cell_list.append((brainArea[i][0], brainArea[i][1]))

		return cell_list

	def pop_cell(self, sessionNr, channelNr_clusterID):
		'''
			This method pops a particular cell to the user
			input:
				sessionNr = the session number that we would like to use
				channelNr_clusterID = the tuple contains both the channel id and the cluster id to select the 
										desired cell
			output:
				cell = a cell object
		'''

		session_name = self.sessions[sessionNr]['session']
		brain_area_file_path = os.path.join(self._construct_data_path(session_name, 'events'), 'brainArea.mat')
		brainArea = loadmat(brain_area_file_path)['brainArea']
		df_brainArea = pd.DataFrame(brainArea)

		brainArea_cell = df_brainArea.loc[(df_brainArea[0] == channelNr_clusterID[0]) & (df_brainArea[1] == channelNr_clusterID[1])]

		cell_path = os.path.join(self._construct_data_path(session_name, 'sorted'), 'A' + str(channelNr_clusterID[0]) + '_cells.mat')
		cell = Cell(cell_path, sessionNr, session_name, *np.asarray(brainArea_cell))
		return cell

	def test(self):
		events = self._get_event_data(25)




class Cell():
	def __init__(self, cell_path, sessionNr, session_name, cell_info):
		self.data_path = cell_path
		self.sessionNr = sessionNr
		self.session_name = session_name
		self.channelNr = cell_info[0]
		self.clusterID = cell_info[1]
		self.oriClusterID = cell_info[2]
		self.brainAreaID = cell_info[3]
		self.spike_timestamps = self._load_cell_data(cell_path)

	def _load_cell_data(self, cell_path):
		'''
			load the raw cell data and capture the spike train timestamps
		'''
		cell_data = loadmat(cell_path)
		channel_spike_timestamps = pd.DataFrame(cell_data['spikes'])
		
		cell_spike_timestamps = np.asarray(channel_spike_timestamps.loc[channel_spike_timestamps[0] == self.clusterID, 2])
		return cell_spike_timestamps





NO_data = NOData()
#print(NO_data.ls_cells(114))
print(NO_data.sessionNrs)
cell_list = NO_data.ls_cells(25)
#print(cell_list)
cell = NO_data.pop_cell(25, (18, 3))
#print(cell.spike_timestamps)
print(cell.spike_timestamps[:])
NO_data.test()
# Test about the indentation size

