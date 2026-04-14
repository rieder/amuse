"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot


#filenames = ["Zeta_M1MSun_dmdt.out", "Zeta_M2MSun_dmdt.out"]
filenames = ["Zeta_M1MSun_dmdt.data"]
ls = ["-", "-", "-", "-.", "--"]
lw = [6, 4, 4, 2, 2]
c = ["b", "k", "g", "b", "r", ]

dmdt = [-10000, -100, -1, -0.01, -0.0001]
zz = 1 # index for the value of zeta
zt = 2# index for the time
zm = 8 # index for the core mass
zd = 12 # index for the value of dmdt
def process_file(filename, mdot):
    t = []
    z = []
    tms = 0
    for line in filename:
        if "mpi" not in line and "Zeta" in line:
            l = line.split()
            if float(l[zd])==mdot:
                #print l
                t.append(float(l[zt]))
                z.append(float(l[zz]))
                if "Main Sequence" in line:
                    tms = t[-1]
    return t, z, tms

def main():
    pyplot.figure()
    pyplot.xlabel("$log_{10}[(t_{end}-t)/t_{MS}]$")
    pyplot.ylabel("$\zeta$")

    pyplot.text(-2.2, 0.8, "giant branch")
    pyplot.text(-0.1, 0.8, "main sequence")
    pyplot.xlim(0.5, -4.)
    pyplot.ylim(-0.6, 1.2)
    fii = 0
    ti=0
    for dmi in dmdt:
        for fi in filenames:
            f = open(fi)
            for line in f:
                t, z, tms = process_file(f, dmi)
                #print t, z, tms
                for i in range(len(t)):
                    if z[i]<0:
                        ti = i
                        break
                for i in range(len(t)):
                    t[i] = numpy.log10((t[-1]-t[i])/tms)
                #            pyplot.scatter(t, z)
                pyplot.plot(t[ti:], z[ti:], ls=ls[fii], lw=lw[fii])
            f.close()
        fii += 1
    pyplot.show()
    #pyplot.savefig("fig_stellar_massloss_response")


if __name__ in ('__main__', '__plot__'):
    main()


