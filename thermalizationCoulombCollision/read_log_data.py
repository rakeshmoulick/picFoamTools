import numpy as np 
import matplotlib.pyplot as plt 
import collision 


#  Creator: Dr. Rakesh Moulick, CPP-IPR

import numpy as np 
import matplotlib.pyplot as plt 
import collision
# Read the data from source
with open('./log.txt','r') as file_handle:
    raw_data = file_handle.readlines()

result_Time = []
result_e    = []
result_Ar   = []
# ---------------------------------------------------
data_time = [row.find("Time") for row in raw_data] 
for i in range (len(data_time)):
    if data_time[i]==0:        
        Time_data = raw_data[i].strip().split(' ')        
        result_Time.append(float(Time_data[2]))
       
# ---------------------------------------------------
data_e = [row.find("    [e]") for row in raw_data]
for i in range (len(data_e)):
    if data_e[i]==0:       
        e_data = raw_data[i].split(' ')        
        result_e.append(float(e_data[37]))            
       
# ---------------------------------------------------
data_Ar = [row.find("    [Ar+]") for row in raw_data]
for i in range (len(data_Ar)):
    if data_Ar[i]==0:       
        Ar_data = raw_data[i].split(' ')
        result_Ar.append(float(Ar_data[37])) 

# print(result_Time)
# print(result_e)
# print(result_Ar)                

result_Time = [i/1e-9 for i in result_Time]

plt.plot(result_Time, result_e, label='$e$')
plt.plot(result_Time, result_Ar, label='$Ar+$')
plt.axhline(15, color='grey', linestyle='--')
plt.xlabel('Time[ns]', fontsize = 15)
plt.ylabel('Temperature[eV]', fontsize = 15)
plt.xticks(fontsize = 15) 
plt.yticks(fontsize = 15) 
plt.legend(loc=0, fontsize=15)
plt.grid(True)
plt.savefig('thermalize.jpg',dpi=300)
plt.show()