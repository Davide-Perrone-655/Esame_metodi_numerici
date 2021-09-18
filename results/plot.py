import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#plt.rcParams.update({'font.size': 18})
'''
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
'''

datas = ['10', '20', '25', '30', '35', '40','45' ,'50','55', '60']

graph = {'ene' : {'fig' : [], 'ax' : [] }, 'mag': {'fig' : [], 'ax' : [] }, 'chi': {'fig' : [], 'ax' : [] }, 'cap': {'fig' : [], 'ax' : [] }}
data = {'ene' : {}, 'mag': {}, 'chi': {}, 'cap': {}}

titles = {'ene' : 'energia', 'mag': 'magnetizzazione', 'chi': 'suscettività magnetica' , 'cap': 'capacità termica'}

for keys in graph:
    graph[keys]['fig'] = plt.figure(figsize = (8, 6), dpi=92)
    graph[keys]['ax'] = plt.subplot(1, 1, 1)
    

print(os.listdir())

for i in datas:
    b, data['ene']['val'], data['ene']['err'],data['mag']['val'], data['mag']['err'], data['chi']['val'], data['chi']['err'], data['cap']['val'], data['cap']['err'] = np.loadtxt('s_res' + i + '.dat', unpack = True)
    for keys in graph:
        graph[keys]['ax'].errorbar(b, data[keys]['val'], data[keys]['err'], marker = '.', linestyle = '', label = i)
        graph[keys]['ax'].set_title('Grafico '+ titles[keys] + ' al variare di L')
        graph[keys]['ax'].grid(True)
        graph[keys]['ax'].legend()
        graph[keys]['ax'].title.set_fontsize(18)
        graph[keys]['ax'].xaxis.label.set_fontsize(18)
        graph[keys]['ax'].yaxis.label.set_fontsize(18)

        graph[keys]['ax'].yaxis.set_minor_locator(AutoMinorLocator(5))
        graph[keys]['ax'].xaxis.set_minor_locator(AutoMinorLocator(5))
        #graph[keys]['ax'].xaxis.grid(True, which='minor')
        
        graph[keys]['fig'].savefig(keys+ '.pdf')




plt.show()

    