# All analyzing tools

class MaxCorrelationClassifier:
    def __init__(self):
        pass
    def fit(self, X, y):
        # senaty check:
        if len(y) != X.shape[1]:
            print('The number of records of X must equal the number of labels')

