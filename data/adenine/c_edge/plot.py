import sys
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.interpolate import spline
from pylab import *
from scipy.interpolate import interp1d

xshift = - 0.5 #=291.24-291.50

data0 = np.genfromtxt('exp.csv')
x0 = data0[:,0] 
y0 = data0[:,1]

data1 = np.genfromtxt('xps.txt')
stcksx1 = data1[:,0] 
stcksy1 = data1[:,1]

data1 = np.genfromtxt('xps.lor')
x1 = data1[:,0] 
y1 = data1[:,1]

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

plt.title('')
plt.xlim(297,288) #Adenine
plt.xlabel('Binding energy (eV)')
yip = 3.9
plt.ylabel('                                            Intensity (arb. units)')
plt.ylim(0.0,yip)#
#ax1.plot(x1+xshift, y1, '#CC0000', label='CVS-EOM-CCSD ( %.2f eV)' %xshift, linewidth=2)#, linestyle='--')
ax1.plot(x1+xshift, y1, 'firebrick', label='fc-CVS-EOM-CCSD ( %.2f eV)' %xshift, linewidth=2)#, linestyle='--')
ax2.plot(x0, y0*4, 'k', label='Experiment', linewidth=2)#, linestyle='--')
n=len(stcksx1)
for i in range(n):
    ax1.plot([stcksx1[i]+xshift,stcksx1[i]+xshift],[0,stcksy1[i]],'firebrick',linewidth=1.0)
#    ax1.plot([stcksx1[i]+xshift,stcksx1[i]+xshift],[0,stcksy1[i]],'#CC0000',linewidth=1.0)
ax1.legend(loc='best',fontsize='small')
ax2.legend(loc='best',fontsize='small')

f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.savefig('adenine_c_xps.pdf')
plt.show()
