import os 
import numpy as np

title = '#b	ene_v	ene_r	mag_v	mag_r	chi_v	chi_r	cap_v	cap_r'

for dfile in os.listdir():
    if dfile.endswith('.dat') and dfile.startswith('res'):
        beta, e_v, e_r, m_v, m_r, c_v, c_r, cp_v, cp_r = np.loadtxt(dfile, unpack = True)
        A_z = sorted(zip(beta, e_v, e_r, m_v, m_r, c_v, c_r, cp_v, cp_r), key=lambda pair: pair[0])
        file1 = open('s_res' + dfile.lstrip('res'), 'w')
        print(title, file=file1)
        temp_1 = [0]*9
        for i in A_z:
            w = ['{:.6f}'.format(j) for j in i]
            if list(i)[0] == temp_1[0]:
                print(list(i))
                print(temp_1)
                continue
            print('\t'.join(w), file = file1 )
            temp_1 = i
        file1.close()
        #temp, dtemp, btemp = [x for x, y, z in sorted(zip(temp, dtemp, betas['40']), key=lambda pair: pair[2])], [y for x, y, z in sorted(zip(temp, dtemp, betas['40']), key=lambda pair: pair[2])], [z for x, y, z in sorted(zip(temp, dtemp, betas['40']), key=lambda pair: pair[2])]
        

