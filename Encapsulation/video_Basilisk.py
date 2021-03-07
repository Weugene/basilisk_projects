# %matplotlib inline
# %matplotlib inline
import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

def gettingFacets(filename, tracer):
    print('Getting facets values')
    if tracer == 1:
        exe = ["./getFacet1", filename]
    else:
        exe = ["./getFacet2", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2,z2)))
                    segs.append(((-r1, z1),(-r2,z2)))
                    skip = True
    return segs


def gettingfield(filename):
    print('Field values')
    exe = ["./getData", filename, str(zmin), str(zmax), str(0), str(rmax),
           str(nx), str(ny)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Rtemp = []
    Ztemp = []
    f1temp = []
    f2temp = []
    veltemp = []
    Omegatemp = []

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            f1temp.append(float(temp3[2]))
            f2temp.append(float(temp3[3]))
            veltemp.append(float(temp3[4]))
            Omegatemp.append(float(temp3[5]))
    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    f1 = np.asarray(f1temp)
    f2 = np.asarray(f2temp)
    vel = np.asarray(veltemp)
    Omega = np.asarray(Omegatemp)
    R.resize((nx+1, ny+1))
    Z.resize((nx+1, ny+1))
    f1.resize((nx+1, ny+1))
    f2.resize((nx+1, ny+1))
    vel.resize((nx+1, ny+1))
    Omega.resize((nx+1, ny+1))
    R.transpose()
    Z.transpose()
    f1.transpose()
    f2.transpose()
    vel.transpose()
    Omega.transpose()

    return R, Z, f1, f2, vel, Omega

# ----------------------------------------------------------------------------------------------------------------------


nGFS = 527

folder = 'VideoBas'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

LEVEL = 10
rmin = -2.32
rmax = 2.32
zmin = -2.16
zmax = 2.16
nx = 2**(LEVEL)
ny = 2**(LEVEL)
lw = 4

if not os.path.isdir(folder):
    os.makedirs(folder)
for ti in range(nGFS):
    t = 0.01 * ti
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%4.4d.png" %(folder, ti)
    if not os.path.exists(place):
        print("%s File not found!" % place)
    else:
        if os.path.exists(name):
            print("%s Image present!" % name)
        else:

            segs1 = gettingFacets(place,1)
            segs2 = gettingFacets(place,2)
            if (len(segs1) == 0 or len(segs2) == 0):
                print("Problem in the available file %s" % place)
            else:
                print("Doing %s" % place)
                R, Z, f1, f2, vel, Omega = gettingfield(place)

                ## Part to plot
                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)

                rc('axes', linewidth=2)
                plt.xticks(fontsize=35)
                plt.yticks(fontsize=35)
                ax.set_xlabel(r'$\mathcal{R}$', fontsize=50)
                ax.set_ylabel(r'$\mathcal{Z}$', fontsize=50)
                ax.set_title('t = %5.4f' % t, fontsize=30)

                facets1=np.array([117,112,179])/255.
                facets2=np.array([217,95,2])/255.

                if ti == 0:
                    Vmax = 1
                else:
                    Vmax = vel.max()

                cntrl1 = ax.imshow(vel, interpolation='bilinear', cmap="hot_r",
                                   origin='lower',
                                   extent=[0, rmax, zmin, zmax],
                                   vmax = Vmax, vmin = 0)
                ## |V|
                cntrl2 = ax.imshow(Omega, interpolation='bilinear', cmap="coolwarm",
                                   origin='lower',
                                   extent=[0, -rmax, zmin, zmax],
                                   vmax = abs(Omega).max(), vmin = -abs(Omega).max())


                ## Drawing Facets
                line_segments1 = LineCollection(segs1, linewidths=lw, colors=facets1, linestyle='solid')
                ax.add_collection(line_segments1)
                line_segments2 = LineCollection(segs2, linewidths=lw, colors=facets2, linestyle='solid')
                ax.add_collection(line_segments2)

                ax.plot([0, 0], [zmin, zmax],'--',color='grey',linewidth=3)

                ax.set_aspect('equal')
                ax.set_xlim(-2.1, 2.1)
                ax.set_ylim(2.1, -2.1)
                ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                ax.yaxis.set_major_locator(plt.MaxNLocator(5))

                l, b, w, h = ax.get_position().bounds

                cb1 = fig.add_axes([l+0.5*w+0.025*w, 1-0.04, 0.5*w-0.1*w, 0.03])
                c1 = plt.colorbar(cntrl1,cax=cb1,orientation='horizontal')
                c1.set_label(r'$\|V\|$',
                             fontsize=30,labelpad=5)
                c1.ax.tick_params(labelsize=20)
                c1.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
                cb2 = fig.add_axes([l+0.025*w, 1-0.04, 0.5*w-0.1*w, 0.03])
                c2 = plt.colorbar(cntrl2,cax=cb2,orientation='horizontal')

                c2.ax.tick_params(labelsize=20)
                c2.set_label(r'$\omega$',fontsize=30,labelpad=5)
                c2.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))

                plt.savefig(name, bbox_inches="tight")
                plt.close()
