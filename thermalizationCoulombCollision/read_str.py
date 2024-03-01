import numpy as np 
import matplotlib.pyplot as plt 

# create a list of names and marks
list1 = ['sravan', 98, 'harsha', 'jyothika', 
         'deepika', 78, 90, 'ramya']

# Get the indicies of the list items containing a string
for i in list1:
    if (type(i) is str):
        j = list1.index(i)
        #print(j, list1[j])


# Get the list index based on a particular string item
list2 = ['sravan kumar', 'rakesh moulick', 'uttam kumar', 'sravan joshi', 'malay aditya']

# indx will list all the instances where the word sravan is present
# it will have a value 0 if the word is present, otherwise -1. 
indx = [row.find("sravan") for row in list2]

for i in range (len(indx)):
    if indx[i]==0:
        print(i, list2[i])
   
