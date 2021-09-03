import numpy as np
import matplotlib.pyplot as plt
from numpy.random import Generator, PCG64
import typing as tp


seed = 10

def setseed(a):
	'''Imposta il seed ad un valore scelto. Default: orario attuale'''

	global seed
	seed=a
	print('seed: %d' %seed)



rg = Generator(PCG64(seed))



def bootstrap(osservabile, vec: tp.List[float], boot_cycle: int = 20) -> float:
    '''Bootstrap function with increasing binning'''
    bootk = []
    len_vec = len(vec)
    #Starts with bin > 1, to reach the asymptotic behaviour faster
    bin_vec = 1 + len_vec//10000

    #Max number of iterations: 7
    while(bin_vec <= 1 + len_vec/10):
        Obs=[]
        #Resampling
        for _ in range(boot_cycle):
            temp = []
            for _ in range(int(len_vec/bin_vec)):
                i = int(rg.random()*len_vec)
                temp.extend(vec[i:min(i+bin_vec,len_vec)])
                
            Obs.append(osservabile(temp))
        bootk.append(np.std(Obs))

        #Increasing binning exponentially
        bin_vec*=2
    return bootk



def media_abs(vec: tp.List[float], valass: bool = False) -> float:
    '''Calculates mean value of a vector, with absolute value of elements if required'''
    if valass:
        return np.mean([abs(i) for i in vec])
    else:
        return np.mean(vec)


def varianza_abs(vec: tp.List[float], valass: bool = False) -> float:
    '''Calculates standard deviation of a vector, with absolute value of elements if valass=True'''
    if valass:
        return np.var([abs(i) for i in vec])
    else:
        return np.var(vec)


def binder(vec: tp.List[float]) -> float:
    '''Calculates binder cumulant'''
    m2=0
    m4=0
    for i in vec:
        m2+=i**2/len(vec)
        m4+=i**4/len(vec)

    if m2 == 0:
        print('Error: division by zero encountered.\nThe current simulation has not enough steps (all zero values picked in binder or bootstrap).\nReturned default expected value (high temperature and L) for binder cumulant: 3')
        return 3

    return m4/m2**2




		
		
	
#nuova funzione: data blocking		
def blocking(A, N=1, jump=0, R_print=False, *args, **kwargs):
	'''funzione di data blocking: calcola l'errore con blocchi di larghezza sempre maggiore, senza ripescare
	A: vettore di dati
	N: dimensione del blocco come 2**N
	jump: cose da saltare in caso di mancate termalizzazioni
	R_print: true per stampare il grafico con le varie iterazioni
	
	'''
	New=(len(A)-jump)
	Div=2**N
	
	if (New<=Div):							#controllo per evitare di tagliare in modo esagerato
		d=int(np.log2(New/1000))
		print('N troppo grande, provare %d' %d)
		return
	
	L1=New - (New%Div)
	B=np.zeros(L1)
	B[:]=A[ jump : (jump + L1)]
	
	mean1=np.mean(B)
	
	if (R_print==True):
		Saved_data=np.zeros(N)
		
		for i in range(N):
			var1=0
			inc=2**(i+1)
			L2=L1//inc					#riduco il mio array, lo preparo per il bootstrap
			
			for j in range(L2):				#medio i valori: definisco il mio vettore con bin
				temp1=(sum(B[j:j+inc]))/inc
				var1=var1 + (((temp1-mean1)**2)/L2)
				
			Saved_data[i]=np.sqrt(var1/(L2-1))
			
		X=np.arange(0,N)
		plt.errorbar(X, Saved_data, markersize=5, linestyle='', marker='.')
		plt.title("Errore da data blocking")
		plt.yscale('log')
		plt.ylabel("Errore")
		plt.xlabel('Lunghezza bin, 2**N')
		plt.grid(color = 'silver')
		plt.show()
		
		res=Saved_data[N-1]
		
	else:
		var1=0
		inc=2**N
		L2=L1//inc					#riduco il mio array, lo preparo per il bootstrap
		for j in range(L2):				#medio i valori: definisco il mio vettore con bin
			temp1=(sum(B[j:j+inc]))/inc
			var1=var1 + (((temp1-mean1)**2)/L2)
				
		res=np.sqrt(var1/(L2-1))
		
	return res
