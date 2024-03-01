#  Creator: Dr. Rakesh Moulick, CPP-IPR

import numpy as np 
import matplotlib.pyplot as plt 
import collision
# Read the data from source
with open('./log.txt','r') as file_handle:
    raw_data = file_handle.readlines()

a = 78 # last line of Time minus one
b = a + 18 #87 

result_Time = []
result_e = []
result_Ar = []

while a <= 317549: # last line of Time minus one    
    key_data = raw_data[a:b]
    key_data = np.array(key_data)
    key_data = key_data.T
        
    # --------------------------------------------
    # Parse Data
    key_data = [row.strip() for row in key_data]      

    data_Time = key_data[0].split(' ') 
    data_e  = key_data[len(key_data)-3].split(' ')
    data_Ar = key_data[len(key_data)-2].split(' ')   
    
    print(a, data_e)
    print(len(key_data))
   
    # ----------------------------------------
    result_Time.append(float(data_Time[2]))    
    result_e.append(float(data_e[30]))    
    result_Ar.append(float(data_Ar[30]))
    a += 21
    b += 21
# ----------------------------------------
# print(result_Time)
# print(result_Ne)
# print(result_Ar)

result_Time = [row for row in result_Time]

plt.plot(result_Time, result_e, color = 'r', label='$T_{Ne}$')   
plt.plot(result_Time, result_Ar, color = 'b', label='$T_{Ar}$') 

# t, y = collision.coll()
# plt.plot(t/1E-6,y[:,0], color = 'r', linestyle='--', label='$T_{Ne (theory)}$')
# plt.plot(t/1E-6,y[:,1], color = 'b', linestyle='--', label='$T_{Ar (theory)}$')

plt.xlabel('Time [ns]')
plt.ylabel('Temperature [eV]')
# plt.xlim(0, 2)
# plt.ylim(500, 1000)
plt.legend(loc=0)
plt.grid(True)
plt.savefig('thermalize.jpg', dpi=300)
plt.show()
