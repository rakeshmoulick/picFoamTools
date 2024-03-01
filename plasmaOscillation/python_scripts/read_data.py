import numpy as np 
import matplotlib.pyplot as plt 

ke = np.loadtxt("ke_data.txt", unpack=True)
# print(len(ke))

ke_sum = 0

ke_sum = np.sum(ke)
print(ke_sum)
