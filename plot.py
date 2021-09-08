import numpy as np
import matplotlib.pyplot as plt

datas = ['20', '30', '40']

graph = {'ene' : {'fig' : [], 'ax' : [] }, 'mag': {'fig' : [], 'ax' : [] }, 'chi': {'fig' : [], 'ax' : [] }, 'cap': {'fig' : [], 'ax' : [] }}
data = {'ene' : {}, 'mag': {}, 'chi': {}, 'cap': {}}


for keys in graph:
    graph[keys]['fig'] = plt.figure(figsize = (8, 6), dpi=92)
    graph[keys]['ax'] = plt.subplot(1, 1, 1)

for i in datas:
    b, data['ene']['val'], data['ene']['err'],data['mag']['val'], data['mag']['err'], data['chi']['val'], data['chi']['err'], data['cap']['val'], data['cap']['err'] = np.loadtxt('res' + i + '.dat', unpack = True)
    for keys in graph:
        graph[keys]['ax'].errorbar(b, data[keys]['val'], data[keys]['err'], marker = '.', markersize = 0.3, linestyle = '')
        graph[keys]['ax'].set_title(keys)
plt.show()
    