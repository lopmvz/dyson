#!/usr/env python
import math
from sys import argv, stdout
import sys
import numpy as np
from pylab import *

def lorentz(omega, exci, sigma,forza):
    grecopi=3.14159265358979323846264338327950288
    return (sigma*forza)/(grecopi*((omega-exci)**2+sigma**2))

def read_file( fd, st, sn ):
    nstate = 0
    state = {}
    energy = {}
    left_norm = {}
    right_norm = {}
    left_dyson = {}
    right_dyson = {}
    lineit = iter(fd)
    for line in lineit:
        if 'Basis set in general basis input format' in line:
            curline = next(lineit)
            curline = next(lineit)
            curline = next(lineit)
            curline = next(lineit)
            natoms = 1 
            file = open("basis", "w")
            file.write("%i    0" %natoms)
            file.write("\n")
            while (not "end" in curline):
                if '****' in curline:
                    curline = next(lineit)
                    if (not 'end' in curline):
                        natoms += 1
                        curline = next(lineit)
                        file.write("\n")
                        file.write("%i    0" %natoms)
                        file.write("\n")
                else:
                    file.write(curline)
                    curline = next(lineit)
            file.close()
        if 'EOM-EE-CCSD state %s -- EOM-IP-CCSD state' %st in line:
        #if 'EOM-EE-CCSD state 1/A\' -- EOM-IP-CCSD state' in line:
        #if 'EOM-EE-CCSD state 1/A\" -- EOM-IP-CCSD state' in line:
            nstate +=1
            state[nstate-1] = (line.split()[4])
            curline = next(lineit)
            energy[nstate-1] = float(curline.split()[3])
            curline = next(lineit)
            curline = next(lineit)
            left_norm[nstate-1] = float(curline.split()[5])
            while (not "Decomposition over AOs for the left alpha Dyson orbital" in curline):
                curline = next(lineit)
            left_dyson[nstate-1] = {}
            i = 0
            curline = next(lineit)
            while (not "*****" in curline):
                left_dyson[nstate-1][i] = float(curline.split()[0])
                i += 1
                curline = next(lineit)
            while (not "Right Dyson orbital norm" in curline):
                curline = next(lineit)
            right_norm[nstate-1] = float(curline.split()[5])
            while (not "Decomposition over AOs for the right alpha Dyson orbital" in curline):
                curline = next(lineit)
            right_dyson[nstate-1] = {}
            i = 0
            curline = next(lineit)
            while (not "*****" in curline):
                right_dyson[nstate-1][i] = float(curline.split()[0])
                i += 1
                curline = next(lineit)
         
        if 'MOLDEN-FORMATTED INPUT FILE FOLLOWS' in line:
            curline = next(lineit)
            file = open("geo", "w")
            while (curline != "\n"):
                file.write(curline)
                curline = next(lineit)
            file.close()

    fgeo = open("geo", "r")
    geo = fgeo.read()
    fbasis = open("basis", "r")
    basis = fbasis.read()
    file1 = open("xps.txt", "w")
    file3 = open("xps.lyp", "w")
    for j in range(len(state)):
        file1.write("{0:.4f}".format((energy[j])))
        file3.write("{0:.2f}".format((energy[j])))
        file1.write("    ")
        file3.write(" & ")
        file1.write("{0:.8f}".format(math.sqrt(left_norm[j]*right_norm[j])))
        file3.write("{0:.6f}".format(math.sqrt(left_norm[j]*right_norm[j])))
        file1.write("\n")
        file3.write(" \\\\ \n")
        nstate = j+1
        file2 = open("left_dyson_%i" %nstate, "w") 
        file2.write(geo)
        file2.write("[GTO]")
        file2.write("\n")
        file2.write(basis)
        file2.write("\n")
        file2.write("[MO]\n")
        file2.write("Sym=X\n")
        file2.write("Ene=")
        file2.write(str(left_norm[j]))
        file2.write("\n")
        file2.write("Spin=Alpha\n")
        file2.write("Occup=1\n")
        k = 0
        for k in range(len(left_dyson[j])):
            file2.write(str(k+1))
            file2.write("    ")
            file2.write(str(left_dyson[j][k]))
            file2.write("\n")
        file2.write("[5D]")

def convolute():
    data = np.genfromtxt("xps.txt")
    exci_SD1 = data[:,0] 
    fCCSD1 = data[:,1]
    file = open("xps.lor","w")
    exci_auSD1=exci_SD1*(1.0/27.2114)
    gamma = 0.4 #eV
    sigma = gamma/2
    n=len(exci_auSD1)
    #for i in range(n):
    #    fCCSD1[i] = fCCSD1[i]#*3/(2*exci_auSD1[i])
    step = 0.01
    omega = np.arange(min(exci_SD1)-5,max(exci_SD1)+5,step) 
    speCCSD1 = []
    for i in range(0, len(omega), 1):
           speCCSD1.append(0.0)
    for i in range(0, len(exci_SD1), 1):
        for n in range(0, len(omega), 1):
            speCCSD1[n]=speCCSD1[n]+lorentz(omega[n],exci_SD1[i],sigma,fCCSD1[i])
    for i in range(len(omega)):
        file.write("{0:.4f}".format(omega[i]))
        file.write("    ")
        file.write("{0:.6f}".format(speCCSD1[i]))
        file.write("\n")

def plot():
    xshift = - 0. #=291.24-291.50

    #data0 = np.genfromtxt('adenine_c_exp.csv')
    #x0 = data0[:,0] 
    #y0 = data0[:,1]
    
    data1 = np.genfromtxt('xps.txt')
    stcksx1 = data1[:,0] 
    stcksy1 = data1[:,1]
    
    data1 = np.genfromtxt('xps.lor')
    x1 = data1[:,0] 
    y1 = data1[:,1]
    
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    
    plt.title('')
#    plt.xlim(560,530) 
#    plt.xlim(306,285) 
    plt.xlabel('Excitation energy (eV)')
    yip = 2
    plt.ylabel('Intensity (arb. units)')
    #plt.ylim(0.0,yip)#
    plt.plot(x1+xshift, y1, 'red', label='CVS-EOM-CCSD ( %.2f eV)' %xshift, linewidth=2)#, linestyle='--')
    n=len(stcksx1)
    for i in range(n):
        plt.plot([stcksx1[i]+xshift,stcksx1[i]+xshift],[0,stcksy1[i]],'r-',linewidth=1.0)
    plt.legend(loc='best',fontsize='small')
    
    #f.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    
    plt.savefig('uracil_xps_%s.pdf' %sn)
    plt.show()

if __name__ == '__main__':
    with open(argv[1],'r') as fd:
        #st = "1/A\""
        #sn = "s1"
        st = "1/A\'"
        sn = "s2"
        read_file(fd,st,sn)
        convolute()
        plot()
