import numpy as np


data = np.loadtxt('/d/scratch/ASTR5160/week13/line.data', unpack=True)
print(len(data))

# KFH The covariance matrix should be 10x10 because we are looking
# KFH at the relatedness of the different bins, of which there are 10,
# KFH not the relatedness of individual randomly generated datapoints in a bin
cov = np.cov(data)

var = [np.var(col, ddof=1) for col in data]
diag  = np.diag(cov)

print("Variance: ", var)
print("Diagonal of covariance matrix:", diag)

cor = np.corrcoef(data)

# KFH Find the location and value of the minimum and maximum non-1 
# KFH correlations between columns of data
mins = [[np.where(cor == i.min()), min(i)] for i in cor]
maxes = [[np.where(cor == max(i[i<.9999])), max(i[i<.9999])] for i in cor]

# KFH The most correlated columns are column 5 and 8 (index 4 and 7)
print("Minimum (loc, value)", mins)
print("Maximum (loc, value)", maxes)
