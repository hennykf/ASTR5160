import numpy as np
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt
from sklearn import neighbors

# KFH Most of this code directly adapted from Adam Myers' example code

iris = load_iris()
colnames = "sepal_length sepal_width petal_length petal_width" 

# KFH Plot sepal length vs width
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for i in range(3):
    target_class = iris.target == i
    ax.scatter(iris.data[target_class, 0], iris.data[target_class, 1],
               s=90, label=iris.target_names[i])

    ax.tick_params(labelsize=14)
    ax.set_xlabel("Sepal length (cm)", size=14)
    ax.set_ylabel("Sepal width (cm)", size=14)
    ax.legend(prop={'size': 14})
fig.savefig('/d/www/kianah/public_html/week11/iris_scatter.png')

# KFH Use k nearest neighbors algorithm with distances only to nearest neighbor
knn = neighbors.KNeighborsClassifier(n_neighbors=1)
knn.fit(iris.data, iris.target)
mock_data = [5, 4, 1, 0]
print(knn.predict([mock_data]), iris.target_names[knn.predict([mock_data])])
mock_data = [6, 3, 4, 1]
print(knn.predict([mock_data]), iris.target_names[knn.predict([mock_data])])

# KFH Map out sepal length/width space
n=100000
mock_data = []

for i in range(2):
    print("working on column: {}".format(colnames.split()[i]))
    col_min = np.min(iris.data[..., i])
    col_max = np.max(iris.data[..., i])
    # ADM generate random points in the space corresponding to the
    # ADM iris measurement of interest.
    mock_meas = np.random.random(n)*(col_max - col_min) + col_min
    mock_data.append(mock_meas)

# KFH reshape array into row-col format instead of list
mock_data = np.reshape(mock_data, (2, n)).T

# KFH Classify new irises based on real world data
knn = neighbors.KNeighborsClassifier(n_neighbors=1)
knn.fit(iris.data[..., :2], iris.target)
mock_target_class = knn.predict(mock_data)

# KFH Print results of classification
for i in range(10):
    print(mock_data[i], mock_target_class[i], iris.target_names[mock_target_class[i]])



# KFH Plot the sepal_length, sepal_width space of classifier
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for i in range(3):
    target_class = mock_target_class == i
    ax.scatter(mock_data[target_class, 0], mock_data[target_class, 1],
               s=10, label=iris.target_names[i])
    ax.tick_params(labelsize=14)
    ax.set_xlabel("Sepal length (cm)", size=14)
    ax.set_ylabel("Sepal width (cm)", size=14)
    ax.legend(prop={'size': 14})

fig.savefig('/d/www/kianah/public_html/week11/iris_zones.png')

print("Percent of virginica iris: ", len(mock_target_class[mock_target_class==2])/100000)



