import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import csv

filename = 'ol_data.csv'
read = csv.reader( open(filename,'r',newline=''), delimiter = ',')

col = next(read)
n_col = len(col)
N = sum(1 for row in read) 

r_roots = np.array([0.155051,0.644949,1.0])
t = np.zeros(N)

print(r_roots)

k      = np.zeros(N) 
cp     = np.zeros(N)
lt1    = np.zeros(N)
vb1    = np.zeros(N)
d1     = np.zeros(N)
b1     = np.zeros(N)
lt2    = np.zeros(N)
vb2    = np.zeros(N)
d2     = np.zeros(N)
b2     = np.zeros(N)
m1_F   = np.zeros(N)
m1_1   = np.zeros(N)
m1_NT  = np.zeros(N)
x1_b_1 = np.zeros(N)
x1_b_2 = np.zeros(N)
x1_c_1 = np.zeros(N)
x1_c_2 = np.zeros(N) 

read = csv.reader( open(filename,'r',newline=''), delimiter = ',')
print(next(read))

i = 0
for r in read:
    k[i]         = float(r[0] )
    cp[i]        = int(r[1] )
    lt1[i]       = float(r[2] )
    vb1[i]       = float(r[3] )
    d1[i]        = float(r[4] )
    b1[i]        = float(r[5] )
    lt2[i]       = float(r[6] )
    vb2[i]       = float(r[7] )
    d2[i]        = float(r[8] )
    b2[i]        = float(r[9] ) 
    m1_F[i]      = float(r[10])
    m1_1[i]      = float(r[11])
    m1_NT[i]     = float(r[12])
    x1_b_1[i]    = float(r[13])
    x1_b_2[i]    = float(r[14])
    x1_c_1[i]    = float(r[15])
    x1_c_2[i]    = float(r[16])
    i = i+1 


for i in range(0,N):
    c = int(cp[i]-1)
    t[i] = k[i] + r_roots[c] 

#plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 16

x1B_b_ref = 1.281271e-02*np.ones(N)
x1A_b_ref = 3.219625e-01*np.ones(N)
m1_F_ref = 5.833503e-01*np.ones(N)
m1_1_ref = 5.347829e-01*np.ones(N)

fig2 = plt.figure(2)
plt.subplot().tick_params(direction='in')
plt.plot(t,x1_b_1,label='mol fraction benzene')
plt.plot(t,x1_b_2,label='mol fraction toluene')
plt.plot(t,x1B_b_ref)
plt.plot(t,x1A_b_ref)
plt.legend()
plt.xlim((0,100))
plt.ylabel('mole fractions')
plt.xlabel('time (min)')
fig2.savefig('boiler_comp.png')

fig3 = plt.figure(3)
plt.subplot().tick_params(direction='in')
plt.plot(t,m1_F_ref)
plt.plot(t,m1_1_ref)
plt.plot(t,m1_1,label='reboiler')
plt.plot(t,m1_F,label='feed tray')
plt.xlim((0,100))
plt.ylabel('holdup (kmol)')
plt.xlabel('time (min)')
fig3.savefig('holdups.png')

#fig4 = plt.figure(4)
#plt.subplot().tick_params(direction='in')
#plt.plot(t,m1_f,label='M1[f]')
#plt.plot(t,m2_f,label='M2[f]')
#plt.plot(t,m1_ntm1,label='M1[NT-1]')
#plt.plot(t,m2_ntm1,label='M2[NT-1]')
#plt.legend()
#plt.ylabel('holdup (kmol)')
#plt.xlabel('time (min)')
#plt.xlim((0,25))
#fig4.savefig('m.png')
#
#fig5 = plt.figure(5) 
#plt.subplot().tick_params(direction='in')
#plt.plot(t,f)
#plt.xlim((0,25))
#plt.ylabel('flow rate (kmol/min)')
#plt.xlabel('time (min)')
#fig5.savefig('f.png')
#
#fig6, ax1 = plt.subplots()
#ax1.plot(t,m2_ntm1,label='M2[NT-1]')
#ax1.set_xlabel('time (min)')
#ax1.set_ylabel('holdup (kmol)')
#ax1.tick_params(direction='in')
#ax1.set_xlim((0,25))
#h1, l1 = ax1.get_legend_handles_labels()
#
#ax2 = ax1.twinx()
#ax2.set_ylabel('flow rate (kmol/min)')
#ax2.plot(t,vb2,label='VB2',color='g')
#ax2.tick_params(direction='in')
#h2, l2 = ax2.get_legend_handles_labels()
##ax2.legend(h1+h2,l1+l2)
#fig6.savefig('state_full_noleg.png',bbox_inches='tight')
#
#fig7, ax1 = plt.subplots()
#ax1.plot(t,err_x,label='x')
#ax1.set_xlabel('time (min)')
#ax1.set_ylabel('x error')
#ax1.tick_params(direction='in')
#ax1.set_xlim((0,25))
#h1, l1 = ax1.get_legend_handles_labels()
#
#ax2 = ax1.twinx()
#ax2.set_ylabel('M error')
#ax2.plot(t,err_m,label='M',color='g')
#ax2.tick_params(direction='in')
#h2, l2 = ax2.get_legend_handles_labels()
#ax2.legend(h1+h2,l1+l2)
#fig7.savefig('err_full.png',bbox_inches='tight')
