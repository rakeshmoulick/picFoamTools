#  Creator: Dr. Rakesh Moulick, CPP-IPR

import numpy as np 

# Read the data from source
with open('./log.txt','r') as file_handle:
    raw_data = file_handle.readlines()

a = 69
b = 87 
key_data = raw_data[a:b]
key_data = np.array(key_data)
key_data = key_data.T

for i in range(len(key_data)):
    print('%s'%(key_data[i]))
print("=======================\n") 
# --------------------------------------------
# Parse Data
key_data = [row.strip() for row in key_data]
for i in range(len(key_data)):
    print('%s'%(key_data[i]))    
# ========================Separate Method Trial (Non Working Yet) =============================================
# for i in range(len(key_data)):    
#     print(i, key_data[i])
# key_data = [row.strip().replace(' ',' ') for row in key_data]
# print(key_data)

# for row in key_data:
#     if (row =='[Ne]                             = 995.4518453 K == 0.08578128489 eV'):
#         data = row
#         #print(row)
#         elif(row == '')

# print(data)

# data_Ne = data.split(' ')
# print(data)
# result_Ne = []
# result_Ne.append(float(data_Ne[30]))
# print(result_Ne)
# =====================================================================
# for i in range(len(key_data)):
#     key_data[i] = [row.strip(" ") for row in key_data[i]]
   
# =====================================================================
print("========================")
print(key_data[0])
print(key_data[15])
print(key_data[16])

print("=======================")  
data_Time = key_data[0].split(' ') 
data_Ne = key_data[15].split(' ')
data_Ar = key_data[16].split(' ')
print(data_Time)
print(data_Ne)
print(data_Ar)
print("=======================")
result_Time = []
result_Time.append(float(data_Time[2]))

result_Ne = []
result_Ne.append(float(data_Ne[30]))

result_Ar = []
result_Ar.append(float(data_Ar[30]))

# for i in range(len(result)):
print(result_Time)
print(result_Ne)
print(result_Ar)