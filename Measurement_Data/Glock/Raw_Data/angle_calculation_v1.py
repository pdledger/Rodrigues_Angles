################################################################################
#This program calculates the complex angle of two complex vectors.
#
#Input data from the eigenvector output file of MPT measurement software.
#This verison of code - input data file is added in the script.
#
#Eigenvectors are of the form: (a + bj)x + (c + dj)y + (e + fj)z : this refers
#to the main data columns combining 1 and 4, 2 and 5, 3 and 6   (excluding freq
#data and balnk padding)
#
#Input file reading strips out blank lines and blank spaces between data.  Each
#frequency has three eigenvectors.
#
#Angles calculated by dot product divided by norm product of between two eigenvectors
#resulting in three complex angles (for angles between Eigenvectors 1,2  1,3 and 2,3.
#
#Output data written to Output.csv
#
#                                               JL Davidson, last update: 9/1/24
#
################################################################################

# Importing
import sys
from time import time
import numpy as np
import subprocess
import cmath
import os
from numpy.linalg import norm
import csv
import array

#open file with eigenvectors - file has header, blank lines and spaces (nulls) in csv format
#file open reads in all data and saves to D1_list with header, blank lines and nulls removed 
with open('Test_eigenvectors.csv') as f:
  reader = csv.reader(f)
  next(reader) #dummy read of header

  D1_list = [] #create data list
  for row in reader:
      if any(row): #pick up non-blank row
         if len(row) != 1: #only append to list the non-blank rows
            this_row = row #read current row
            new_this_row = ' '.join(this_row).split() #remove nulls in data
            D1_list.append(new_this_row) #add data with removed nulls to list 



#convert list data to array - data now comprises of 7 column-wide data of:
# freq  EigVect(Re)1, x, y, z components, EigVect(Im)1, x, y, z components             
# freq  EigVect(Re)2, x, y, z components EigVect(Im)2, x, y, z components    
# freq  EigVect(Re)3, x, y, z components, EigVect(Im)3, x, y, z components    
# .... then the next 3 rows at the next frequency step
D2_array = np.array(D1_list, dtype=float)

#separate eigenvectors and frequencies for later processing
Eig_vect = D2_array[:,1:7] #Eigenvectors (3 Re components and 3 Im components)
F1 = D2_array[:,0] #frequencies with repeats
F1_length = len(F1) #array length
Freq = F1[0:F1_length:3] #access every 3rd element to remove repeated frequencies

#print to screen frequencies and the number of frequency points
print(' ')
print(Freq)
print("\nLength of frequency array...\n",len(Freq))
print(' ')

#loop thru' Eigenvector array and combine the real and imaginary parts
row = len(Eig_vect) #number of rows
Evec_comp1 = [] #initialise array
for i in range(0,row): #per  row
   EvRe_1 = Eig_vect[i][0] #1st component in Re Eigenvectors 1,2,3 (all frequencies)
   EvRe_2 = Eig_vect[i][1] #2nd component in Re Eigenvectors 1,2,3 (all frequencies)
   EvRe_3 = Eig_vect[i][2] #3rd component in Re Eigenvectors 1,2,3 (all frequencies)
   EvIm_1 = Eig_vect[i][3] #1st component in Im Eigenvectors 1,2,3 (all frequencies)
   EvIm_2 = Eig_vect[i][4] #2nd component in Im Eigenvectors 1,2,3 (all frequencies)
   EvIm_3 = Eig_vect[i][5] #3rd component in Im Eigenvectors 1,2,3 (all frequencies)

   #combine Re and Im parts for the components of Eigenvectors (all frequencies) 
   Evect1 = [EvRe_1 + 1j*EvIm_1] 
   Evect2 = [EvRe_2 + 1j*EvIm_2]
   Evect3 = [EvRe_3 + 1j*EvIm_3]

   #combine Evect1, Evect2 and Evect3 into one array
   this_Evec_comp1 = Evect1 + Evect2 + Evect3
   Evec_comp1.append(this_Evec_comp1) #append to output array
   
#rearrange Evec_comp1 to group eigenvectors row by row for a single frequency
#output data to new array Evec_compl_Reshape
k = int((row / 3)) #takes into account each frequency is repeated for 3 rows for eignevectors associated at that frequency
cnt = 0 #set counter at first row
Evec_comp1_Reshape = [] #initialise array
for p in range(0,k):

   EigVect_1 = Evec_comp1[cnt][0:5] #pull out first eigenvector data
   EigVect_2 = Evec_comp1[cnt + 1][0:5] #pull out second eigenvector data
   EigVect_3 = Evec_comp1[cnt + 2][0:5] #pull out third eigenvector data

   #write to new array and reformat data - each row is for a single frequency
   #of 3 eigenvectors (each with x, y and z real and imaginary components)
   EigVectGrp = EigVect_1 + EigVect_2 + EigVect_3 
   Evec_comp1_Reshape.append(EigVectGrp) #append to output array 

   cnt = cnt + 3 #update counter to next frequency group
   

#print to screen reshaped complex eigenvector data
print(' ')
print(Evec_comp1_Reshape) 
print("\nLength of eigenvector array...\n",len(Evec_comp1_Reshape))
print(' ')


#calculation of angles between eigenvectors
#angles calculated by dot product divided by norm product
#three output angles refer to:
# theta1: angle between EigVector1 and EigVector2
# theta2: angle between EigVector1 and EigVector3
# theta3: angle between EigVector2 and EigVector3
#each outputted angle is complex

print (len(Evec_comp1_Reshape[0]))


AllThetaData = [] #initialise array
for q in range(0,len(Evec_comp1_Reshape)):


    #define each eigenvector
    v1 = Evec_comp1_Reshape[q][0:3]
    v2 = Evec_comp1_Reshape[q][3:6]
    v3 = Evec_comp1_Reshape[q][6:9]

    #create numpy arrays
    arr1 = np.array(v1)
    arr2 = np.array(v2)
    arr3 = np.array(v3)

    #display the current vectors for calculation 
    print(' ')
    print("\neigenvector 1 ...\n",arr1)
    print("\neigenvector 2 ...\n",arr2)
    print("\neigenvector 3 ...\n",arr3)

    #calculate the vector dot product 
    vdot12 = np.vdot(arr1,arr2) #for vectors 1 and 2
    vdot13 = np.vdot(arr1,arr3) #for vectors 1 and 3
    vdot23 = np.vdot(arr2,arr3) #for vectors 2 and 3

    #display result to screen
    print("\nDot Product - eigenvectors 1 and 2 ...\n",vdot12)
    print("\nDot Product - eigenvectors 1 and 3 ...\n",vdot13)
    print("\nDot Product - eigenvectors 2 and 3 ...\n",vdot23)

    #calculate norm products of complex vectors (for vectors 1 and 2)
    norm_prod12 = norm(v1)*norm(v2)

    #calculate norm products of complex vectors (for vectors 1 and 3)
    norm_prod13 = norm(v1)*norm(v3)

    #calculate norm products of complex vectors (for vectors 2 and 3)
    norm_prod23 = norm(v2)*norm(v3)

    #display result to screen
    print("\nNorm Product - eigenvectors 1 and 2 ...\n",norm_prod12)
    print("\nNorm Product - eigenvectors 1 and 3 ...\n",norm_prod13)
    print("\nNorm Product - eigenvectors 2 and 3 ...\n",norm_prod23)

    #calculate complex angles  
    theta1 = cmath.acos(vdot12 / norm_prod12) #angle between eigenvectors 1 and 2
    theta2 = cmath.acos(vdot13 / norm_prod13) #angle between eigenvectors 1 and 3
    theta3 = cmath.acos(vdot23 / norm_prod23) #angle between eigenvectors 2 and 3

    #print angles to screen
    print("\nCalculated Angle, theta1 ...\n",theta1)
    print("\nCalculated Angle, theta2 ...\n",theta2)
    print("\nCalculated Angle, theta3 ...\n",theta3)

    #Group together all theta data and write to new array
    ThetaGrp = theta1, theta2, theta3 
    AllThetaData.append(ThetaGrp) #append to output array 

#print output angle data to screen
print(' ')
print(AllThetaData)     
print(len(AllThetaData))

#combine Freq data with calculated angle data
#for writing to screen
#OutputData = np.column_stack((Freq, AllThetaData))
#print(OutputData)



#write data to csv file (no freq. data)
header = ['Theta1', 'Theta2', 'Theta3']
with open('Output.csv', 'w', newline='') as ff:
    writer = csv.writer(ff, delimiter=',')
    writer.writerow(header)
    writer.writerows(AllThetaData)

print("\nData written to Output.csv\n")
    

