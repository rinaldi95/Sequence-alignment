import numpy as np
from numba import njit, cuda,vectorize,prange
import time
import glob
from multiprocessing import Pool
import multiprocessing
import math as m
import threading

'''Fill the matrix except for last line and column, the default value for the diagonal
is =mismatch, if there's a match the value become =match.
Choose the max value from the diagonal+match/mismatch, vertical+indel and orizontal+indel.'''
@njit
def matrix_generator(vertical,nw,s1,s2):
  
    match=20
    transizione=-10
    trasversione=-30
    open_gap=(-80) #Insertion or deletion
    length_gap=(-10)
   
    for i in range(len(nw[0])):
        nw[0][i]=i*length_gap
    for i in range (len(vertical)):
        vertical[i]+=open_gap+(i)*length_gap
   
    for i in range (1,len(s1)+1):
        
        #if i%500==0:
        #    print(i)
        idisparity=(i-1)%2
        nw[idisparity][0]=(i-1)*length_gap
        
        horizontal=((i-1)*length_gap)+open_gap

        for j in range (1, len(nw[0])):
            jminus1=j-1
            
           
            add=trasversione
            if s1[i-1]==s2[j-1]:
                add=match
            if ((s1[i-1]=='A' and s2[j-1]=='G')or 
                (s1[i-1]=='G' and s2[j-1]=='A')or
                (s1[i-1]=='T' and s2[j-1]=='C')or
                (s1[i-1]=='C' and s2[j-1]=='T')):
                add=transizione
                
            vertical[jminus1]+=length_gap
            horizontal+=length_gap

            max_list=[vertical[jminus1], horizontal, nw[idisparity][jminus1]]
            nw[i%2][j]=add+max(max_list)

            if (nw[idisparity][jminus1]+open_gap>=horizontal):
                horizontal=nw[idisparity][jminus1]+open_gap
  
            if (nw[idisparity][jminus1]+open_gap>=vertical[jminus1]):
                vertical[jminus1]=nw[idisparity][jminus1]+open_gap
           
    return nw[i%2][j]

#@njit(parallel=True)
def smith_waterman(distance_matrix,list_genoma,n_process,cores,list_files):
    for line in range (len(list_genoma)):
        for column in range (line+1,len(list_genoma)):
            if ((line*len(list_genoma)+column-int((line+1)*(line+2)/2))%cores)==n_process:
                

                s1=np.array(list_genoma[line])
                s2=np.array(list_genoma[column])
                
                
                nw=np.zeros((2,len(s2)+1),dtype=int)
                
                
                vertical=np.zeros(len(nw[0]),dtype=int)
             
                t1=time.time()
                
                score=matrix_generator(vertical,nw,s1,s2)
                print('process n°%i: genoma %i vs %i completed, time:%f'%(n_process,line,column,time.time()-t1))
                
                distance_matrix[line*len(list_genoma)+column]=score
                distance_matrix[column*len(list_genoma)+line]=score
                del nw
                
               
    print('process %i terminated'%n_process)


if __name__=='__main__':
    list_files=glob.glob('p2\\*.gb')#
    list_files=sorted(list_files)
    print(list_files)

    list_genoma=[]
    for i in range( len(list_files)):
        genoma=[]
        f = open(list_files[i], "r")                            #open the file. remeber to specify the directory if you are not inside the same directory of the file.
        x='a'                                       #inizializzation of x is required if you want to compare: (x[0:6]!='ORIGIN') and (x!='')

        if f.read(1)=='>':                              #all FASTA file starts with ">"
            print('file format: FASTA')
    #        g=open(sys.argv[2],"w+")                       #open destination file
            f.readline()                                #send empty a line
            for x in f.read():                          #for every single char after the first line
                if x in ('T','G','C','A'):            #what is 'N'? there are other possibilities that i'm ignoring?
                    genoma.append(x)
    #                print(x,end='')                        #print in terminal. you can comment this line. end='' avoid the newline
    #               g.write(x)                      #print in file
    #        g.close()                              #close the file

        elif f.read(4)=='OCUS':                             #all Gencode file starts with "LOCUS". i have just read a character so only "OCUS" remain
            print('file format: Gencode')
            while (x[0:6]!='ORIGIN') and (x!=''):                   #the nucleic acid sequesnce start with "ORIGIN". x=='' when you arrive at end of file
                x=f.readline()                          #now x is a list so i can write (x[0:6]!='ORIGIN')
                x.strip                             #deleting some spaces

            if x=='':
                print('i haven\'t found the word \"ORIGIN\"')

            else:
    #            g = open(sys.argv[2],"w+")                 #same as FASTA format
                while(f.readline(2)!='//'):                 #the nucleic acid sequence end with "//"
                    for x in f.readline():
                        if x in ('t','g','c','a'):
                            genoma.append(x)
    #                       print(x.upper(),end='')         #I want the same output for FASTA and Gencode. so i make an uppercase character
    #                       g.write(x.upper())
    #           g.close()
        else:
            print('il file non è nè in formato FASTA nè in formato Gencode')    #if you don't find ">" or "LOCUS"

        f.close
        list_genoma.append(genoma)

    list_genoma=np.array(list_genoma,dtype=object)
    
    #distance_matrix=np.zeros((len(list_files),len(list_files)),dtype=int)

    distance_matrix=multiprocessing.Array('i',len(list_files)**2,lock=False)
    cores=int(multiprocessing.cpu_count())
    
    
    print('utilizzerò %i core'%cores)
    process=[multiprocessing.Process(target=smith_waterman, args=(distance_matrix,list_genoma,i,cores,list_files)) for i in range (cores)]
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
    f=open("distance_matrix_genoma_2010308010.ods","w+")
    for i in range (len(list_files)):
        f.write(','+list_files[i].replace('.gb',''))
    f.write('\n')
    for i in range (len(list_files)):
        f.write((list_files[i].replace('.gb',''))+',')
        for j in range (len(list_files)):
            f.write("%i,"%distance_matrix[i*len(list_files)+j])
        f.write("\n")
    f.close()
    del distance_matrix





