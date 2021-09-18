import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
from copy import deepcopy
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def f1(x, a, b, c):
    return a*(x**2) + b*x +c 

def p_max(a, b, c):
    return -b/(2*a), (c - (b**2/(4*a)))

def f2(L, a, b, c):
    return a + b*(L**c)

beta_crit = 0.4407

datas = ['10', '20', '25', '30', '35', '40','45' ,'50','55', '60']

graph = {'ene' : {'fig' : [], 'ax' : [] }, 'mag': {'fig' : [], 'ax' : [] }, 'chi': {'fig' : [], 'ax' : [] }, 'cap': {'fig' : [], 'ax' : [] }}
data_2 = {'ene' : {}, 'mag': {}, 'chi': {}, 'cap': {}}

betas = {}
res = {}

for i in datas:
    res[i] = deepcopy(data_2)
    betas[i], res[i]['ene']['val'], res[i]['ene']['err'],res[i]['mag']['val'], res[i]['mag']['err'], res[i]['chi']['val'], res[i]['chi']['err'], res[i]['cap']['val'], res[i]['cap']['err'] = np.loadtxt('s_res' + i + '.dat', unpack = True)
    




#stima brutale
chi_m = [np.max(res[i]['chi']['val']) for i in betas]
dchi_m = [np.max(res[i]['chi']['err']) for i in betas]

b_c = [betas[i][list(res[i]['chi']['val']).index(np.max(res[i]['chi']['val']))] for i in betas]
#
chi_m1 =[]
b_c1 = []
dchi_m1 = []
db_c1 = []


####################################à
#   fit picchi suscettività
rg = 4

for n in datas:
    
    max_1 = list(res[n]['chi']['val']).index(np.max(res[n]['chi']['val']))
    par_max = {'x': [], 'y': []}


    for i in range(rg):
        plt.errorbar(betas[n], res[n]['chi']['val'], res[n]['chi']['err'], linestyle = '', marker = '.')
        temp1 = res[n]['chi']['val'][max_1 - (i+3) : max_1 + (i+2)]
        dtemp1 = res[n]['chi']['err'][max_1 - (i+3) : max_1 + (i+2)]
        betas1 = betas[n][max_1 - (i+3) : max_1 + (i+2)]
        plt.errorbar(betas1, temp1, dtemp1, linestyle = '', marker = '.')
        popt, pcovm = curve_fit(f1, betas1, temp1, sigma = dtemp1 )
        print(popt)
        t_1, t_2 = p_max(*popt)
        par_max['x'].append(t_1)
        par_max['y'].append(t_2)
        x_axis = np.linspace(betas[n][max(max_1 - (i+2), 0)], betas[n][min(max_1 + (i+2), len(betas[n]))], 1000)
        #plt.plot(x_axis, f1(x_axis, *popt))
        #plt.show()
    #plt.errorbar(b_c, chi_m, dchi_m, linestyle='', marker = '.' )
    print(par_max)
    chi_m1.append(np.mean(par_max['y']))
    b_c1.append(np.mean(par_max['x']))
    #dchi_m1.append(np.sqrt(np.std(par_max['x'])**2 + np.std(par_max['y'])**2 ))
    dchi_m1.append(np.std(par_max['y']))
    db_c1.append(np.std(par_max['x']))
    plt.show()
    plt.close()
    
###################
#   plot figure varie 
fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
ax.errorbar(b_c, chi_m, dchi_m, linestyle = '', marker = '.', color = 'red')
ax.errorbar(b_c1, chi_m1, dchi_m1, linestyle = '', marker = '.', color = 'blue')
ax.grid(True, color = 'silver')
#plt.show()
plt.close()




####################
#   stima gamma 
L_axis = [int(i) for i in datas]


fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
ax.errorbar(L_axis, chi_m1, dchi_m1, linestyle = '', marker = '.', color = 'green')
ax.grid(True, color = 'silver')

in_val = [0, 0.1, 7/4]
x_axis = np.linspace(L_axis[0], L_axis[-1], 1000)

#ax.plot(x_axis, f2(x_axis, *in_val))


popt, pcovm = curve_fit(f2, L_axis, chi_m1, p0=in_val, sigma = dchi_m1)
print(popt)
print(np.sqrt(np.diag(pcovm)))
gamma = popt[2]
dgamma = np.sqrt(np.diag(pcovm))[2]
print('esponente critico gamma: {} +- {}'.format(popt[2], np.sqrt(np.diag(pcovm))[2]))
ax.plot(x_axis, f2(x_axis, *popt))
#plt.show()
plt.close()

#########################
#stima nu

fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
cut = 1
ax.errorbar(L_axis[cut:], b_c1[cut:], db_c1[cut:], linestyle = '', marker = '.', color = 'orange')
ax.grid(True, color = 'silver')

in_val = [0.44, -0.1, -1]
x_axis = np.linspace(L_axis[0], L_axis[-1], 1000)

#ax.plot(x_axis, f2(x_axis, *in_val))
#plt.show()


popt, pcovm = curve_fit(f2, L_axis[cut:], b_c1[cut:], p0=in_val, sigma = db_c1[cut:])
print(popt)
print(np.sqrt(np.diag(pcovm)))
nu = popt[2]
dnu = np.sqrt(np.diag(pcovm))[2]
print('esponente critico nu: {} +- {}'.format(popt[2], np.sqrt(np.diag(pcovm))[2]))
ax.plot(x_axis, f2(x_axis, *popt))
#plt.show()
plt.close()


#########################
#plot fss

fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)

valid_markers = ([item[0] for item in mpl.markers.MarkerStyle.markers.items() if item[1] is not 'nothing' and not item[1].startswith('tick') and not item[1].startswith('caret')])
v_mark = ['.', '^', '*', 'X']
k = len(valid_markers) -1


for i in datas:
    beta_plot = [(j - beta_crit)*(int(i)**(-nu)) for j in betas[i]]
    chi_plot = [j/(int(i)**gamma) for j in res[i]['chi']['val']]
    dchi_plot = [j/(int(i)**gamma) for j in res[i]['chi']['err']]
    ax.errorbar(beta_plot, chi_plot, dchi_plot, linestyle = '', marker = valid_markers[k], label = i)
    k-=1

ax.set_title('Finite size scaling suscettività')
ax.grid(True, color = 'silver')

ax.legend()
ax.title.set_fontsize(18)
ax.xaxis.label.set_fontsize(18)
ax.yaxis.label.set_fontsize(18)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

fig.savefig('fss_chi.pdf')
plt.show()

#print(betas)
#print(res)

#for keys in graph:
#    graph[keys]['fig'] = plt.figure(figsize = (8, 6), dpi=92)
#    graph[keys]['ax'] = plt.subplot(1, 1, 1)
#
#    graph[keys]['ax'].errorbar(b, data[keys]['val'], data[keys]['err'], marker = '.', markersize = 5, linestyle = '')
#    graph[keys]['ax'].set_title(keys)
#    graph[keys]['ax'].grid(True)