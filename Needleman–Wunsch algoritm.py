import numpy as np
from numba import njit, cuda,vectorize,prange
import time

'''Fill the matrix except for last line and column, the default value for the diagonal
is =mismatch, if there's a match the value become =match.
Choose the max value from the diagonal+match/mismatch, vertical+indel and orizontal+indel.'''

@njit
def matrix_generator(nw,match,mismatch,indel,start,step):
    
    for i in range (start, min(len(nw)-1,start+step)):
        
        for j in range (1, len(nw[0])-1):
            add=mismatch                                        
            if s1[i-1]==s2[j-1]:                                    
                add=match                                   
            nw[i][j]=max(nw[i][j-1]+indel, nw[i-1][j]+indel, nw[i-1][j-1]+add)  

    return nw


s2="acactgatcg" 
s1="acactg"   

nw=np.zeros((len(s1)+1,len(s2)+1),dtype=int) 

#Score assignment 

match=1                                                 
mismatch=0
indel=-1    #Insertion or deletion

t1=time.time()
step=100
for i in range (0,len(nw)-1,step):
    print (i)
    nw=matrix_generator(nw,match,mismatch,indel,i,step)
print('tempo impiegato=%f'%(time.time()-t1))

#Fill the last line and column avoiding to add indel     
j=len(s2)
for i in range(1,len(nw)-1):
    add=mismatch
    if s1[i-1]==s2[j-1]:
        add=match
    nw[i][j]=max(nw[i][j-1]+indel, nw[i-1][j], nw[i-1][j-1]+add)    #Add indel only horizontally
i=len(s1)
for j in range(1,len(nw[0])-1):
    add=mismatch
    if s1[i-1]==s2[j-1]:
        add=match
    nw[i][j]=max(nw[i][j-1], nw[i-1][j]+indel, nw[i-1][j-1]+add)    #Add indel only vertically

#Fill the last position avoiding to add indel at all
j=len(s2)
add=mismatch
if s1[i-1]==s2[j-1]:
    add=match
nw[i][j]=max(nw[i][j-1], nw[i-1][j], nw[i-1][j-1]+add)

#Print the matrix in a file
f=open("matrix.txt","w+")
for i in range (len(nw)):
    for j in range (len(nw[0])):
        f.write("%i "%nw[i][j]) 
    f.write("\n")
f.close()

#Here the alligned sequences will be wrote
a1=[]   
a2=[]

# Trace back to origin the path that produce the max score.
# Again there are different rules for the edges of the matrix.
# Since the first control is along the diagonal, if there are more possibilities,
# the priority will be given in this order: diagonal, vertical, horizontal

#Start from the end
i=len(s1)
j=len(s2)

#First position

#Calculating the diagonal score
add=mismatch
if s1[i-1]==s2[j-1]:
    add=match

#If the score is due to a diagonal movement, move diagonally and append to  
#the arrays the next nucleic acid from both the sequences
if (nw[i][j]==nw[i-1][j-1]+add):    
    i=i-1
    j=j-1
    a1.append(s1[i])
    a2.append(s2[j])

#If the score is due to a vertical or diagonal movement, move diagonally and
#append to an array the next nucleic acid for the sequences along wich you move
#and a "-" for the sequence along wich you didn't move 
elif (nw[i][j]==nw[i-1][j]):
    i=i-1
    a1.append(s1[i])
    a2.append("-")
elif (nw[i][j]==nw[i][j-1]):
    j=j-1
    a1.append("-")
    a2.append(s2[j])

#Repeat the same operation for all the matrix.
#The rules change if you are in the first/last line/column.

#First line
while(i==len(s1))and(j!=len(s2)):
    add=mismatch
    if s1[i-1]==s2[j-1]:
        add=match
    if (nw[i][j]==nw[i-1][j-1]+add):
        i=i-1
        j=j-1
        a1.append(s1[i])
        a2.append(s2[j])
    elif (nw[i][j]==nw[i-1][j]+indel):
        i=i-1
        a1.append(s1[i])
        a2.append("-")
    elif (nw[i][j]==nw[i][j-1]):
        j=j-1
        a1.append("-")
        a2.append(s2[j])

#first column
while(i!=len(s1))and(j==len(s2)):
    add=mismatch
    if s1[i-1]==s2[j-1]:
        add=match
    if (nw[i][j]==nw[i-1][j-1]+add):
        i=i-1
        j=j-1
        a1.append(s1[i])
        a2.append(s2[j])
    elif (nw[i][j]==nw[i-1][j]):
        i=i-1
        a1.append(s1[i])
        a2.append("-")
    elif (nw[i][j]==nw[i][j-1]+indel):
        j=j-1
        a1.append("-")
        a2.append(s2[j])
    
#Centre of the matrix
while(i!=0)and(j!=0):
    add=mismatch
    if s1[i-1]==s2[j-1]:
        add=match
    if (nw[i][j]==nw[i-1][j-1]+add):
        i=i-1
        j=j-1
        a1.append(s1[i])
        a2.append(s2[j])
    elif (nw[i][j]==nw[i-1][j]+indel):
        i=i-1
        a1.append(s1[i])
        a2.append("-")
    elif (nw[i][j]==nw[i][j-1]+indel):
        j=j-1
        a1.append("-")
        a2.append(s2[j])

#Last line. Don't need controls. Just go to nw[0][0] putting all "-"
while(i!=0):
    i=i-1
    a1.append(s1[i])
    a2.append("-")  

#Last column.
while(j!=0):    
    j=j-1
    a1.append("-")
    a2.append(s2[j])
    
#The alligned sequences are in reverse. Need to be inverted.            
a1.reverse()
a2.reverse()

#Finally print the final output
f=open("sars mers.txt","w+")
f.write("score:"+str(nw[len(nw)-1][len(nw[0])-1])+'\n')

for i in range (len(a1)):
    
    f.write("%s %s"%(a1[i],a2[i]))
    if a1[i]!=a2[i]:
        if a1[i]=="-" or a2[i]=="-":
            f.write("deletion")
        else: f.write("mismatch")
    
    f.write("\n")

#use this in order to have have lines of 60 char    
# j=0
# k=0
# print (len(a1))        
# for i in range (int(2*len(a1)/120)*120):
    # if i%120<60 and i%120>=0:
        # f.write("%s"%(a1[j]))
        # j+=1
        # print(i)
    # if i%120==59:
        # f.write("\n") 
    # if i%120<120 and i%120>=60:
        # f.write("%s"%(a2[k]))
        # k+=1
    # if i%120==119:
        # f.write("\n\n")        
# while j<len(a1):
    # f.write("%s"%(a1[j]))
    # j+=1
# f.write('\n')
# while k<len(a2):
    # f.write("%s"%(a1[k]))
    # k+=1         
        
f.close()






  
