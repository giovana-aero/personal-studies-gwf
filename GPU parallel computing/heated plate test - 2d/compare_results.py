import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt("results_original.txt",delimiter = " ")
data2 = np.genfromtxt("results.txt",delimiter = " ")

compare = []
for j in range(1,data1.shape[1]-3):
    for i in range(1,data1.shape[0]-2):
        compare.append(abs(data1[j,i]-data2[j,i]))

print(np.mean(compare))