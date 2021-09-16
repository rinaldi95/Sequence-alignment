import numpy as np
from numba import njit, cuda,vectorize,prange
import time
import glob
from multiprocessing import Pool
import multiprocessing
import math as m
import threading
from sklearn import preprocessing

'''Fill the matrix except for last line and column, the default value for the diagonal
is =mismatch, if there's a match the value become =match.
Choose the max value from the diagonal+match/mismatch, vertical+indel and orizontal+indel.'''
@njit
def matrix_generator(start,vertical_pointers,vertical,nw,punt,p,s1,s2):

    #Score assignment
    match=10
    transizione=-15
    trasversione=-30
    open_gap=-100 
    length_gap=-10

    vertical_pointers-=1

#inizialization of the variables that contains the max across the columns
    for i in range (len(nw[0])):
        vertical[i]=length_gap*i+open_gap
        
#inizialization of the edge of the matrix    
    for i in range (len(nw)):
        nw[i][0]=length_gap*i
    for i in range (len(nw[0])):
        nw[0][i]=length_gap*i  
    
    for i in range (start, len(nw)):
#inizialization of the variables that contains the max across this line    
        horizontal=(open_gap)+nw[i-1][0]
        vertical_pointers+=1
#p points to the box wich from we are calculating the score       
        p[1]=1
        
        for j in range (1, len(nw[0])):
            add=trasversione
            if s1[i-1]==s2[j-1]:
                add=match
            if ((s1[i-1]=='A' and s2[j-1]=='G')or 
                (s1[i-1]=='G' and s2[j-1]=='A')or
                (s1[i-1]=='T' and s2[j-1]=='C')or
                (s1[i-1]=='C' and s2[j-1]=='T')):
                add=transizione
                
#I increase the max across this line and column every step                
            vertical[j-1]+=length_gap
            horizontal+=length_gap
            
#and so the distance o f the box from wich this max is calculated
            p[1]-=1
            p[0]=vertical_pointers[j-1]

#I chose the max, find the score, and save the position of the box of origin
            max_list=[vertical[j-1], horizontal, nw[i-1][j-1]]
            nw[i][j]=add+max(max_list)
            punt[i][j]=p[max_list.index(max(max_list))]
           
#I check if there is a new max for this line and column           
            if (nw[i-1][j-1]+open_gap>=horizontal):
                horizontal=nw[i-1][j-1]+open_gap
                p[1]=0
            if (nw[i-1][j-1]+open_gap>=vertical[j-1]):
                vertical[j-1]=nw[i-1][j-1]+open_gap
                vertical_pointers[j-1]=0
    
    return nw,punt

#@njit(parallel=True)
def smith_waterman(distance_matrix,list_genoma,n_process,cores,list_files):
    for line in range (len(list_genoma)):
        for column in range (line+1,len(list_genoma)):
            if ((column+line*len(list_genoma))%cores)==n_process:

                s1=np.array(list_genoma[line])
                s2=np.array(list_genoma[column])

                nw=np.zeros((len(s1)+1,len(s2)+1),dtype=int)
                
                punt=np.zeros((len(s1)+1,len(s2)+1),dtype=int)
                p=np.zeros(3,dtype=int)
                vertical_pointers=np.zeros(len(nw[0]),dtype=int)
                vertical=np.zeros(len(nw[0]),dtype=int)
                
                lenght_nw_col=len(nw)
                lenght_nw_line=len(nw[0])
                t1=time.time()

                
                
                nw,punt=matrix_generator(1,vertical_pointers,vertical,nw,punt,p,s1,s2)
                
                print(nw)
                i=len(nw)-1
                j=len(nw[0])-1
                a1=[]
                a2=[]
                print (punt)
                print(nw[-1][-1])
                while i!=0 and j!=0 :
                    i-=1
                    j-=1
                    #print('i=%f' %i)
                    a1.append(s1[i])
                    a2.append(s2[j])
                    if punt[i+1][j+1]>0:
                        for k in range(punt[i+1][j+1]):
                            a2.append('-')
                            a1.append(s1[i-(k+1)])
                        i-=punt[i+1][j+1]
                    if punt[i+1][j+1]<0:
                        for k in range(-punt[i+1][j+1]):
                            a2.append(s2[j-(k+1)])
                            a1.append('-')
                        j+=punt[i+1][j+1]

                del punt
                
                for k in range(i-1,-1,-1):
                    a1.append(s1[k])
                    a2.append('-')

                for k in range(j-1,-1,-1):
                    a1.append('-')
                    a2.append(s2[k])

                #The alligned sequences are in reverse. Need to be inverted.
                a1.reverse()
                a2.reverse()
                
                        #Finally print the final output
                print('printing the alligned sequencies')
                f=open(list_files[line].replace('.gb','')+(list_files[column].replace('.gb','.txt')).replace('p1\\',''),"w+")
                for i in range (len(a1)):
                    f.write("%s %s\n"%(a1[i],a2[i]))
                f.close()

                distance_matrix[line*len(list_genoma)+column]=nw[-1][-1]
                distance_matrix[column*len(list_genoma)+line]=nw[-1][-1]
                del nw
                
    print('il processo %i è concluso'%n_process)


if __name__=='__main__':
    list_files=glob.glob('*.gb')
    list_files=sorted(list_files)
    print(list_files)

    list_genoma=[]
    for i in range( len(list_files)):
        genoma=[]
        f = open(list_files[i], "r")                #open the file. remeber to specify the directory if you are not inside the same directory of the file.
        x='a'                                       #inizializzation of x is required if you want to compare: (x[0:6]!='ORIGIN') and (x!='')

        if f.read(1)=='>':                          #all FASTA file starts with ">"
            print('file format: FASTA')
    #        g=open(sys.argv[2],"w+")               #open destination file
            f.readline()                            #send empty a line
            for x in f.read():                      #for every single char after the first line
                if x in ('T','G','C','A'):          #what is 'N'? there are other possibilities that i'm ignoring?
                    genoma.append(x)
    #                print(x,end='')                #print in terminal. you can comment this line. end='' avoid the newline
    #               g.write(x)                      #print in file
    #        g.close()                              #close the file

        elif f.read(4)=='OCUS':                     #all Gencode file starts with "LOCUS". i have just read a character so only "OCUS" remain
            print('file format: Gencode')
            while (x[0:6]!='ORIGIN') and (x!=''):   #the nucleic acid sequesnce start with "ORIGIN". x=='' when you arrive at end of file
                x=f.readline()                      #now x is a list so i can write (x[0:6]!='ORIGIN')
                x.strip                             #deleting some spaces

            if x=='':
                print('i haven\'t found the word \"ORIGIN\"')

            else:
    #           g = open(sys.argv[2],"w+")                 #same as FASTA format
                while(f.readline(2)!='//'):                 #the nucleic acid sequence end with "//"
                    for x in f.readline():
                        if x in ('t','g','c','a'):
                            genoma.append(x.upper())
    #                       print(x.upper(),end='')         #I want the same output for FASTA and Gencode. so i make an uppercase character
    #                       g.write(x.upper())
    #           g.close()
        else:
            print('il file non è nè in formato FASTA nè in formato Gencode')    

        f.close
        list_genoma.append(genoma)

    list_genoma=np.array(list_genoma)
    print(len(list_genoma))
    #distance_matrix=np.zeros((len(list_files),len(list_files)),dtype=int)

    distance_matrix=multiprocessing.Array('i',len(list_files)**2,lock=False)
    cores=int(multiprocessing.cpu_count()/2)

    
    print('utilizzerò %i core'%cores)
    process=[multiprocessing.Process(target=smith_waterman, 
        args=(distance_matrix,list_genoma,i,cores,list_files)) for i in range (cores)]
    for i in range (cores):
        print('preparazione processo %i'%i)
        process[i].start()
        print('processo %i avviato'%i)

    for i in range (cores):
        print('attendo il processo %i'%i)
        process[i].join()


    #distance_matrix=np.zeros((len(list_files),len(list_files)),dtype=int)
    #Print the matrix in a file
    print('scrivo la matrice delle distanze su file')
    f=open("distance_matrix.txt","w+")
    for i in range (len(list_files)):
        for j in range (len(list_files)):
            f.write("%i "%distance_matrix[i*len(list_files)+j])
        f.write("\n")
    f.close()
    del distance_matrix


'''
        #Print the matrix in a file
        f=open("scorenew.txt","w+")
        for i in range (len(nw)):
            for j in range (len(nw[0])):
                f.write("%i "%nw[i][j])
            f.write("\n")
        f.close()

        #Print the matrix in a file
        f=open("pointersnew.txt","w+")
        for i in range (len(nw)):
            for j in range (len(nw[0])):
                f.write("%i "%punt[i][j])
            f.write("\n")
        f.close()
        '''
'''
        #Finally print the final output
        f=open(list_files[line]+,"w+")
        for i in range (len(a1)):
            f.write("%s %s\n"%(a1[i],a2[i]))
        f.close()
'''




4