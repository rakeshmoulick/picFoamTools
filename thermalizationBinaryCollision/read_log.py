#  Creator: Dr. Rakesh Moulick, CPP-IPR

import numpy as np 
import matplotlib.pyplot as plt 
import collision
# Read the data from source
with open('./log.txt','r') as file_handle:
    raw_data = file_handle.readlines()

a = 69 # last line of Time minus one
b = a + 18 #87 

result_Time = []
result_Ne = []
result_Ar = []

while a<=4248: # last line of Time minus one
    key_data = raw_data[a:b]
    key_data = np.array(key_data)
    key_data = key_data.T
    # --------------------------------------------
    # Parse Data
    key_data = [row.strip() for row in key_data]      

    data_Time = key_data[0].split(' ') 
    data_Ne = key_data[15].split(' ')
    data_Ar = key_data[16].split(' ')
    # ----------------------------------------
    result_Time.append(float(data_Time[2]))    
    result_Ne.append(float(data_Ne[30]))    
    result_Ar.append(float(data_Ar[30]))
    a += 21
    b += 21
# ----------------------------------------
# print(result_Time)
# print(result_Ne)
# print(result_Ar)

result_Time = [row/1e-6 for row in result_Time]

plt.plot(result_Time, result_Ne, color = 'r', label='$T_{Ne}$')   
plt.plot(result_Time, result_Ar, color = 'b', label='$T_{Ar}$') 
hline = 750
plt.axhline(hline, color='grey', linestyle='--', label='$T_{avg}$')

t, y = collision.coll()
plt.plot(t/1E-6,y[:,0], color = 'r', linestyle='--', label='$T_{Ne (theory)}$')
plt.plot(t/1E-6,y[:,1], color = 'b', linestyle='--', label='$T_{Ar (theory)}$')

plt.xlabel('Time [$\mu$s]')
plt.ylabel('Temperature [K]')
plt.xlim(0, 2)
plt.ylim(500, 1000)
plt.legend(loc=0)
plt.grid(True)
plt.savefig('./thermalize.jpg', dpi=300)
plt.show()
