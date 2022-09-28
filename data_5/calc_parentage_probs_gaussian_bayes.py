import argparse
import numpy as np
from sklearn.naive_bayes import GaussianNB
from sklearn import preprocessing

########
# Args #
########

my_parser = argparse.ArgumentParser(description='Import arguments for the script')

# Add the common arguments
my_parser.add_argument('--training_file', default='data/training_classification.tsv', type=str,
                        help='File with data and correctness of assignments to be used for training')

my_parser.add_argument('--test_file', default=None, help='Data file to assess probability of correct parentage assignments')

args = my_parser.parse_args()

########
# Main #
########

training_data_array = np.genfromtxt(args.training_file, delimiter='\t')

y_train=training_data_array[1:,-1]
X_train=training_data_array[1:,0:2]

print(X_train)
sc = preprocessing.MinMaxScaler().fit(X_train)
X_train_scaled = sc.transform(X_train)


X_test = np.genfromtxt(args.test_file, delimiter=',')
print(X_test)
#X_test_single =  X_test.reshape(1, -1)
X_test_scaled = sc.transform(X_test)


gnb = GaussianNB()

y_pred = gnb.fit(X_train_scaled, y_train).predict(X_test_scaled)

print(y_pred)

#estimate probability
prob = gnb.predict_proba(X_test_scaled)
print(prob)



prob_list = []
for i, data in enumerate(y_pred):
       print("Chance of correct prediction: ", prob[i,1])
       prob_list.append(prob[i,1])

with open('prob_list.txt', 'w') as f:
    for item in prob_list:
        f.write("%s\n" % item)
