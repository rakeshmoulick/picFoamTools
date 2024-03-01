import numpy as np 
import matplotlib.pyplot as plt 
import collision 


#  Creator: Dr. Rakesh Moulick, CPP-IPR

import numpy as np 
import matplotlib.pyplot as plt 
import collision
# Read the data from source
with open('./test_log.txt','r') as file_handle:
    raw_data = file_handle.readlines()

result_Time = []
result_e    = []
result_Ar   = []
# ---------------------------------------------------
data_time = [row.find("Time") for row in raw_data] 
# print(data_time)
for i in range (len(data_time)):
    if data_time[i]==0:
        # print(i, data_time[i])
        # print(raw_data[i])
        Time_data = raw_data[i].strip().split(' ')
        print(Time_data)
        #print(Time_data[2])        
        result_Time.append(float(Time_data[2]))
        #print("===============================")
# ---------------------------------------------------
data_e = [row.find("    [e]") for row in raw_data]

#print(data_e)
for i in range (len(data_e)):
    if data_e[i]==0:
        # print(i, data_e[i])
        # print(i, raw_data[i])
        e_data = raw_data[i].split(' ')
        #print(e_data)
        #print(float(e_data[34])) 
        #e_data.append(float(raw_data[i][30])) 
        result_e.append(float(e_data[37]))
        
# print(result_e)

# e_data = raw_data[15].split(' ')
# print(e_data)
# print(float(e_data[37]))       
       
# ---------------------------------------------------
data_Ar = [row.find("    [Ar+]") for row in raw_data]
# print(data_Ar)
for i in range (len(data_Ar)):
    if data_Ar[i]==0:
        # print(i, data_Ar[i])
        # print(raw_data[i])
        Ar_data = raw_data[i].split(' ')
        result_Ar.append(float(Ar_data[37])) 

print(result_Time)
print(result_e)
print(result_Ar)                

result_Time = [i/1e-9 for i in result_Time]

plt.plot(result_Time, result_e)
plt.plot(result_Time, result_Ar)
plt.show()