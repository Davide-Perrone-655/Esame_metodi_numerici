import numpy as np 
import os 
from bootstrap import boot_1 as bt
from copy import deepcopy
import matplotlib.pyplot as plt

datas = ['45','55']
#datas = ['10','50', '60']
datas_1 = ['20']
res = {'ene' : {'val' : [], 'err' : [] }, 'mag': {'val' : [], 'err' : [] }, 'chi': {'val' : [], 'err' : [] }, 'cap': {'val' : [], 'err' : [] }}
betas = {}
res_2 = {}
for i in datas:
    res_2[i] = deepcopy(res)
    betas[i] = []

#print(res_2)

for i in datas:
    d = int(i)**2
    os.chdir('./'+ i)
    count = 0
    print(os.listdir())
    for meas in os.listdir():
        betas[i].append(float(meas.strip('.dat').split('_')[2]))

        temp, ene_, mag_ = np.loadtxt(meas, unpack = True)
        ene_1 = ene_[int(len(ene_)/10):]
        mag_1 = mag_[int(len(ene_)/10):]
        res_2[i]['ene']['val'].append(bt.media_abs(ene_1))
        res_2[i]['mag']['val'].append(bt.media_abs(mag_1, True))
        res_2[i]['cap']['val'].append(d*bt.varianza_abs(ene_1))
        res_2[i]['chi']['val'].append(d*bt.varianza_abs(mag_1, True))

        res_2[i]['ene']['err'].append(bt.bootstrap(lambda x: bt.media_abs(x), ene_1)[-1])
        res_2[i]['mag']['err'].append(bt.bootstrap(lambda x: bt.media_abs(x, True), mag_1)[-1])
        res_2[i]['cap']['err'].append(bt.bootstrap(lambda x: d*bt.varianza_abs(x), ene_1)[-1])
        res_2[i]['chi']['err'].append(bt.bootstrap(lambda x: d*bt.varianza_abs(x, True), mag_1)[-1])
        
        #print(betas)
        count +=1
        print('Punto finito, %d', count)
        print(meas)
    os.chdir('..')
    #print(betas[i])
    #print(res_2)
    file1 = open('res' + i + '.dat', 'w')
    file1.write('b\tene_v\tene_r\tmag_v\tmag_r\tchi_v\tchi_r\tcap_v\tcap_r\n')
    for j in range(len(betas[i])):
        file1.write('{:.6f}\t'.format(betas[i][j]))
        for keys in res_2[i]:
            for keys_2 in res_2[i][keys]:
                file1.write('{:.6f}\t'.format(res_2[i][keys][keys_2][j]))

        file1.write('\n')
    file1.close()
    print('Cartella finita')
print(res_2)
print(betas)
print(len(res['mag']['val']))


'''
for keys in res_2[0]:
    graph[keys] 
    fig = plt.figure(figsize = (8, 6), dpi=92)
    ax = plt.subplot(1, 1, 1)
    for i in datas:
    
        plt.errorbar(betas[i], res_2[i][keys]['val'], marker = '.', linestyle = '')
    plt.show()
'''
