import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
from copy import deepcopy
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


def chi2(f, popt, x, y, dy):
    return (((y - f(x, *popt))/dy)**2).sum()

def f1(x, a, b, c):
    return a*(x**2) + b*x +c 

def p_max(a, b, c):
    return -b/(2*a), (c - (b**2/(4*a)))

def f2(L, a, b, c):
    return a + b*(L**c)

def lin(x, a, b):
    return a + b*x

def f3(x, a, b, c):
    return a + b*np.log(x) + 0.5*c*(np.log(x)**2)

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
cap_m = [np.max(res[i]['cap']['val']) for i in betas]
dcap_m = [np.max(res[i]['cap']['err']) for i in betas]

b_c = [betas[i][list(res[i]['cap']['val']).index(np.max(res[i]['cap']['val']))] for i in betas]
#
cap_m1 =[]
b_c1 = []
dcap_m1 = []
db_c1 = []


####################################à
#   fit picchi suscettività
rg = 4

for n in datas:
    
    max_1 = list(res[n]['cap']['val']).index(np.max(res[n]['cap']['val']))
    par_max = {'x': [], 'y': []}


    for i in range(rg):
        plt.errorbar(betas[n], res[n]['cap']['val'], res[n]['cap']['err'], linestyle = '', marker = '.')
        temp1 = res[n]['cap']['val'][max_1 - (i+3) : max_1 + (i+3)]
        dtemp1 = res[n]['cap']['err'][max_1 - (i+3) : max_1 + (i+3)]
        betas1 = betas[n][max_1 - (i+3) : max_1 + (i+3)]
        plt.errorbar(betas1, temp1, dtemp1, linestyle = '', marker = '.')
        popt, pcovm = curve_fit(f1, betas1, temp1, sigma = dtemp1 )
        print(popt)
        t_1, t_2 = p_max(*popt)
        par_max['x'].append(t_1)
        par_max['y'].append(t_2)
        x_axis = np.linspace(betas[n][max(max_1 - (i+2), 0)], betas[n][min(max_1 + (i+2), len(betas[n]))], 1000)
        plt.plot(x_axis, f1(x_axis, *popt))
        #plt.show()
    #plt.errorbar(b_c, cap_m, dcap_m, linestyle='', marker = '.', color = 'orange')
    print(par_max)
    cap_m1.append(np.mean(par_max['y']))
    b_c1.append(np.mean(par_max['x']))
    #dcap_m1.append(np.sqrt(np.std(par_max['x'])**2 + np.std(par_max['y'])**2 ))
    dcap_m1.append(np.std(par_max['y']))
    db_c1.append(np.std(par_max['x']))
    #plt.show()
    plt.close()
    
###################
#   plot figure varie 
fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
ax.errorbar(b_c, cap_m, dcap_m, linestyle = '', marker = '.', color = 'red')
ax.errorbar(b_c1, cap_m1, dcap_m1, linestyle = '', marker = '.', color = 'blue')
ax.grid(True, color = 'silver')
#plt.show()
plt.close()




####################
#   stima alpha
L_axis = [int(i) for i in datas]


fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
ax.errorbar(L_axis, cap_m1, dcap_m1, linestyle = '', marker = '.', color = 'green')
ax.grid(True, color = 'silver')

cut = 1
in_val = [5, 1, 0.5]
x_axis = np.linspace(L_axis[0], L_axis[-1], 1000)

ax.plot(x_axis, f2(x_axis, *in_val))
#plt.show()


popt, pcovm = curve_fit(f2, L_axis[cut:], cap_m1[cut:], p0=in_val, sigma = dcap_m1[cut:])
print(popt)
print(np.sqrt(np.diag(pcovm)))
alpha = popt[2]
dalpha = np.sqrt(np.diag(pcovm))[2]
print('esponente critico alpha: {} +- {}'.format(popt[2], np.sqrt(np.diag(pcovm))[2]))
ax.plot(x_axis, f2(x_axis, *popt))
#plt.show()
plt.close()

#######################################
#   nuova stima alpha

L_axis = [ np.log(int(i)) for i in datas]
lcap = np.log(cap_m1)
dlcap = np.log(dcap_m1)

fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
ax.errorbar(L_axis, lcap, dlcap, linestyle = '', marker = '.', color = 'green')


ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True, color = 'silver')
cut = 3
in_val = [3, 0.1]
x_axis = np.linspace(L_axis[0], L_axis[-1], 1000)

#ax.plot(x_axis, lin(x_axis, *in_val))
#plt.show()


popt, pcovm = curve_fit(lin, L_axis[cut:], lcap[cut:], p0=in_val, sigma = dlcap[cut:])
print(popt)
print(np.sqrt(np.diag(pcovm)))
alpha = popt[1]
dalpha = np.sqrt(np.diag(pcovm))[1]
print('Nuovo esponente critico alpha: {} +- {}'.format(popt[1], np.sqrt(np.diag(pcovm))[1]))
ax.plot(x_axis, lin(x_axis, *popt))
#plt.show()
plt.close()






#######################################
#   ennesima nuova stima alpha (linearizzata)

L_axis = [int(i) for i in datas]


fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
ax.errorbar(L_axis, cap_m1, dcap_m1, linestyle = '', marker = '.', color = 'green')
ax.grid(True, color = 'silver')

cut = 1
in_val = [2, 2, 0.1]
x_axis = np.linspace(L_axis[0], L_axis[-1], 1000)

ax.plot(x_axis, f3(x_axis, *in_val))
#plt.show()


popt, pcovm = curve_fit(f3, L_axis[cut:], cap_m1[cut:], p0=in_val, sigma = dcap_m1[cut:])
print(popt)
print(np.sqrt(np.diag(pcovm)))

chi_2 = (((cap_m1[cut:] - f3(L_axis[cut:], *popt))/dcap_m1[cut:])**2).sum()
print('chiquadro: {}, ndof: {}'.format(chi_2, len(cap_m1)-3))

alpha = popt[2]/popt[1]
dalpha = alpha * ((np.sqrt(np.diag(pcovm))[2]/popt[2]) + (np.sqrt(np.diag(pcovm))[1]/popt[1]))
print('Ennesimo esponente critico alpha: {} +- {}'.format(alpha, dalpha))
ax.plot(x_axis, f3(x_axis, *popt))
plt.show()
plt.close()



#########################
#stima nu

L_axis = [int(i) for i in datas]

fig = plt.figure(figsize = (8, 6), dpi=92)
ax = plt.subplot(1, 1, 1)
cut = 0
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
#v_mark = ['.', '^', '*', 'X']
k = len(valid_markers) -1
alpha = 0 
nu = -1
cap_plot_m = np.max(res['60']['cap']['val'])

for i in datas:

    beta_plot = [(j - beta_crit)*(int(i)**(-nu)) for j in betas[i]]
    cap_plot = [j/(int(i)**alpha) for j in res[i]['cap']['val']] - np.max(res[i]['cap']['val']) + cap_plot_m
    #cap_plot = [j/(int(i)**alpha) for j in res[i]['cap']['val']]
    dcap_plot = [j/(int(i)**alpha) for j in res[i]['cap']['err']]
    ax.errorbar(beta_plot, cap_plot, dcap_plot, linestyle = '', marker = valid_markers[k], label = i)
    k-=1

ax.set_title('Finite size scaling capacità termica')
ax.grid(True, color = 'silver')
ax.legend()
ax.title.set_fontsize(18)
ax.xaxis.label.set_fontsize(18)
ax.yaxis.label.set_fontsize(18)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

fig.savefig('fss_cap.pdf')
plt.show()