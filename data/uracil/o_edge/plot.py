import sys
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.interpolate import spline
from pylab import *
from scipy.interpolate import interp1d

xshift = - 0. #=291.24-291.50

#S0
data0 = np.genfromtxt('s0min/s0/xps.txt')
stcksx0 = data0[:,0] 
stcksy0 = data0[:,1]
data0 = np.genfromtxt('s0min/s0/xps.lor')
x0 = data0[:,0] 
y0 = data0[:,1]

#S1
data1 = np.genfromtxt('s0min/s1/xps.txt')
stcksx1 = data1[:,0] 
stcksy1 = data1[:,1]
data1 = np.genfromtxt('s0min/s1/xps.lor')
x1 = data1[:,0] 
y1 = data1[:,1]
#S1min
'''
data12 = np.genfromtxt('s1min/s1/xps.txt')
stcksx12 = data12[:,0] 
stcksy12 = data12[:,1]
data12 = np.genfromtxt('s1min/s1/xps.lor')
x12 = data12[:,0] 
y12 = data12[:,1]
'''

#S2
data2 = np.genfromtxt('s0min/s2/xps.txt')
stcksx2 = data2[:,0] 
stcksy2 = data2[:,1]
data2 = np.genfromtxt('s0min/s2/xps.lor')
x2 = data2[:,0] 
y2 = data2[:,1]
#S2min
'''
data22 = np.genfromtxt('s2min/s2/xps.txt')
stcksx22 = data22[:,0] 
stcksy22 = data22[:,1]
data22 = np.genfromtxt('s2min/s2/xps.lor')
x22 = data22[:,0] 
y22 = data22[:,1]
'''

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

#f, (plt, ax2) = plt.subplots(2, sharex=True, sharey=True)

plt.title('')
plt.xlim(560,530) 
plt.xlabel('Excitation energy (eV)')
yip = 3
plt.yticks([])
plt.ylabel('Intensity (arb. units)')
plt.ylim(0.0,yip)#
plt.plot(x0+xshift, y0, '#CC0000', label='S$_0$ ', linewidth=2)#, linestyle='--')
plt.plot(x1+xshift, y1, '#000080', label='S$_1$ ', linewidth=2)#, linestyle='--')
plt.plot(x2+xshift, y2, '#006400', label='S$_2$ ', linewidth=2)#, linestyle='--')


n=len(stcksx0)
for i in range(n):
    plt.plot([stcksx0[i]+xshift,stcksx0[i]+xshift],[0,stcksy0[i]],'#CC0000',linewidth=1.0)
n=len(stcksx1)
for i in range(n):
    plt.plot([stcksx1[i]+xshift,stcksx1[i]+xshift],[0,stcksy1[i]],'#000080',linewidth=1.0)
n=len(stcksx2)
for i in range(n):
    plt.plot([stcksx2[i]+xshift,stcksx2[i]+xshift],[0,stcksy2[i]],'#006400',linewidth=1.0)
plt.legend(loc='upper right',fontsize='small', framealpha=1)
#f.subplots_adjust(hspace=0)
#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.savefig('uracil_o_tr_xps.pdf')
plt.show()
