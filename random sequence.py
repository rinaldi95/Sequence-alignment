
from random import seed
from random import randint
# seed random number generator

if __name__=='__main__':
    seed(10)
    # generate some random numbers
    f=open('sequenza random4000 2.RANDOM',"w+")
    lunghezza=4000               
    f.write('LOCUS\nORIGIN\n')               
    for i in range (lunghezza):
        a=randint(0,3)
        if a==0:
            f.write('a')
        if a==1:
            f.write('c')
        if a==2:
            f.write('g')
        if a==3:
            f.write('t')

    f.write('\n//')               
    f.close()
