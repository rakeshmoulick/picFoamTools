
# Source: https://www.youtube.com/watch?v=ZXqnwy0pUeE

import numpy as np 

# Read the data from source
with open('./sales_data.txt','r') as file_handle:
    raw_data = file_handle.readlines()

# print(raw_data)

key_products = raw_data[10:19]
key_products = np.array(key_products)
key_products = key_products.T

for i in range(len(key_products)):
    print('%s'%(key_products[i]))
# ----------------------------------------------------------
# Parse the data : Observe the structure and to extract
# key_products = [row.strip() for row in key_products]

key_products = [row.strip().replace(')','').split('(') for row in key_products]

key_products = [[row[0]] + row[1].split(" ")[:2] for row in key_products]

for i in range(len(key_products)):
    print('%s'%(key_products[i]))

# ------------------------------------------------------------
result = []

# Sort Data
for row in key_products:
    if row[2] == 'thousand':
        result.append((row[0], float(row[1])*1000))
    else:
        result.append((row[0], float(row[1])*1000000))

for i in range(len(result)):
    print(result[i])

# Sorting ...
print("\n")
top5 = sorted(result,key = lambda x: x[1], reverse=True)[0:5]
for i in range(len(top5)):
    print(top5[i])

# Printing out in required format
print("\n")
for i, item in enumerate(top5):
    print("{num}. {item_name}".format(num=i, item_name=item[0]))










# ------------------------------------------------------
# data = ['my','name','is','rakesh']
# data = np.array(data)
# data = data.T
# for i in range(len(data)):
#     print('%s'%(data[i]))
# ------------------------------------------------------
