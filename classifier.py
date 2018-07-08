# The maximum correlation coefficient classifier
import numpy as np
#import pandas as pd


class MaxCorrCoefClf:
    def __init__(self):
        self.templates = np.array([])
        self.unique_labels = np.array([])

    def fit(self, x_train, y):
        """
        Fit the data tot the classifier
        :param x_train: A [n_feature, n_record] numpy matrix
        :param y: A [n_record, 1] label vector
        :return:
        """
        # sanity check:
        if len(y) != x_train.shape[1]:
            print('The number of records of X_train must equal the number of labels')

        self.unique_labels = np.unique(y)
        self.templates = np.empty((x_train.shape[0], 0))
        for i in range(len(self.unique_labels)):
            temp = np.mean(x_train[:, y == self.unique_labels[i]], axis=1, keepdims=True)
            self.templates = np.append(self.templates, temp, axis=1)
            
    def predict(self, x_test):
        """
        Calculate the maximum correlation coefficients for each test record and associate it with the labels
        :param x_test: Test cases
        :return: a list with the classified results.
        """
        if x_test.shape[0] > 1:
            mean_subtracted_templates = self.templates - \
                                         np.tile(np.mean(self.templates, axis=0, keepdims=True),
                                                 [self.templates.shape[0], 1])
            mean_subtracted_x_test = x_test - np.tile(np.mean(x_test, axis=0, keepdims=True), [x_test.shape[0], 1])
            normalization_matrix = np.matmul(np.sqrt(np.matrix(np.diag(np.matmul(mean_subtracted_templates.T,
                                                               mean_subtracted_templates)))).T,
                                             np.sqrt(np.matrix(np.diag(np.matmul(mean_subtracted_x_test.T,
                                                               x_test)))))
            template_corrcoeffs = (np.matmul(mean_subtracted_templates.T, 
                                             mean_subtracted_x_test) / normalization_matrix)

        else:
            template_corrcoeffs = -1 * np.square((np.tile(self.templates, [x_test.shape[0], 1])
                                                  - np.tile(x_test.T, [1, self.templates.shape[1]])))

        max_indexes = self._randmax(template_corrcoeffs.T)
        output_labels = [self.unique_labels[ind] for ind in max_indexes]
        return output_labels

    @staticmethod
    def _randmax(template_corrcoeffs):
        """
        Get random indexes if there are two exact same correlation results for each trial.
        :param template_corrcoeffs: The coefficient matrix that we want to look at.
        :return max_indexes: The classification result
        """
        max_indexes = []
        for i in range(template_corrcoeffs.shape[0]):
            output = template_corrcoeffs[i, :]
            sort_index = np.argsort(output)
            if output[0, sort_index[0, -1]] == output[0, sort_index[0, -2]]:
                indexes = [sort_index[0, -1], sort_index[0, -2]]
                max_indexes.append(indexes[np.random.randint(2)])
            else:
                max_indexes.append(sort_index[0, -1])
        return max_indexes
