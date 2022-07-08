import numpy as np
import pandas as pd
import glob
import os
import re
import copy
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import minimize, rosen, rosen_der, brentq
from scipy.integrate import solve_ivp
from numpy import sin, cos, pi , exp
from scipy.integrate import odeint
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from tsmoothie.smoother import *
from scipy.interpolate import UnivariateSpline
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
import networkx as nx
from scipy import interpolate
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit
from numpy import linalg as LA

SMALL_SIZE = 20
MEDIUM_SIZE = 25
BIGGER_SIZE = 30
iter = 0
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.rcParams['font.size'] = 25
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

# set tight margins
plt.margins(0.015, tight=True)
from scipy.stats import norm

def thickness(Ca):
    return 0.67*Ca**(2./3.)/(1 + 3.35*Ca**(2./3.))

def sort_names(image_files):
    file_names = [os.path.basename(string) for string in image_files]
    times = [(float(re.findall("\d+\.\d+", string)[0]), string) for string in file_names]
    times = sorted(times, key=lambda x: x[0])
    print(f"Time: {times[0][0]} -- {times[-1][0]}")
    image_files = [t[1] for t in times]
    return image_files

def fun_psi(psi, z, s1, s2, l, B):
    X, Sigma = z
    if psi < 1e-4:
        return [0, 2/B]
    Q = sin(psi)/Sigma + s1*s2*X/l**2 - B
    return [sin(psi)/Q, -cos(psi)/Q]

def fun_x(x, z, s1, s2, l, B):
    w, f = z
    return [-(1 + w**2)*(1/f + np.sqrt(1 + w**2)*(-s1*s2*x/l**2 - B)), -w]

def volume(X, Sigma):
    y2 = Sigma**2
    Vd_comp = np.abs(np.pi*integrate.trapezoid(y2, x=X))
    return Vd_comp

def target_fun(X, Sigma, Vd):
    Vd_comp = volume(X, Sigma)
    Phi = np.abs(Vd_comp - Vd)/Vd
    # print(f'Vd_comp={Vd_comp} Vd={Vd} Phi={Phi}')
    return Phi
#
# dX/dpsi=sin(psi)/Q
# dSigma/spsi = -cos(psi)/Q
# Q = sin(psi)/Sigma + s1*s2*X/l**2 - B
# X(0)=Sigma(0)=0

def shape_psi(B_, alpha, s1, s2, l):
    B = B_[0]
    X0 = 0
    Sigma0 = 0
    soln_psi = solve_ivp(fun_psi, (0, alpha), (X0, Sigma0), method='BDF', args=(s1, s2, l, B),
                         min_step=1e-6, max_step = 1e-3, rtol = 1e-15, dense_output=True)
    X_psi = soln_psi.y[0]
    Sigma_psi = soln_psi.y[1]
    return X_psi, Sigma_psi

def compute_shape_psi(B_, alpha, s1, s2, l, Vd):
    X_psi, Sigma_psi = shape_psi(B_, alpha, s1, s2, l)
    return target_fun(X_psi, Sigma_psi, Vd)

def compute_regime(B_, l, props):
    X_psi, Sigma_psi = shape_psi(B_, 0.5*np.pi, props['s1'], props['s2'], l)
    Vd_comp = volume(X_psi, Sigma_psi)
    Ub =
    Ca = props['mu1']*Ub/props['sigma']
    if Vd_comp > props['Vd']:
        return
    Phi = np.abs(Vd_comp - Vd)/Vd
#
#w'=f''=(1 + w**2)*(1/f + np.sqrt(1 + w**2)*(s1*s2*x/l**2 - B))
#f'=w
def shape_full(d_, B, s1, s2, l):
    X_psi, Sigma_psi = shape_psi([B], pi/6, s1, s2, l)
    d = d_[0]
    X0 = -X_psi[-1]
    W0 = (Sigma_psi[-1] - Sigma_psi[-2])/(X_psi[-1] - X_psi[-2])
    f0 = Sigma_psi[-1]
    print(f'X0={X0} f0={f0} W0={W0} d={d}, f()={fun_x(X0, (W0, f0), s1, s2, l, B)}')
    soln_x = solve_ivp(fun_x, (X0, d), (W0, f0), method='BDF', args=(s1, s2, l, B),
                       min_step=1e-6, max_step=1e-3, rtol=1e-15, dense_output=True)
    X_x = -soln_x.t
    W_x = soln_x.y[0]
    Sigma_x = soln_x.y[1]
    X = np.concatenate([X_psi, X_x])
    Sigma = np.concatenate([Sigma_psi, Sigma_x])
    return X_psi, Sigma_psi, X_x, Sigma_x
    # return X, Sigma

def compute_shape_full(d_, B, s1, s2, l, Vd):
    X_psi, Sigma_psi, X_x, Sigma_x = shape_full(d_, B, s1, s2, l)
    X = np.concatenate([X_psi, X_x])
    Sigma = np.concatenate([Sigma_psi, Sigma_x])
    return target_fun(X, Sigma, Vd)

# s1 = -1 for pendant drop
# s2 = 1 for rhod - rhoa > 0
def pendant_drop(props, picScale, mode=None):
    Vd = props['Vd']
    alpha = props['alpha']
    s1 = props['s1']
    s2 = props['s2']
    diam = props['diam']
    a = (3*Vd/(4*np.pi))**(1/3)
    props["a"] = a
    d = a*props['d/diam']
    l = np.sqrt(props['sigma']/((props['rho1'] - props['rho2'])*props['grav']))
    iBo = (l/a)**2
    B = 2/(a*(4/(2 + cos(alpha)**3 - 3*cos(alpha)))**(1/3))
    # B *= 0.9
    print(f'Ba_guess={B*a}, d/diam={d/diam}, a={a}, l/diam={l/diam}, a/Bo={iBo*a}')
    if not mode:
        # res = minimize(compute_shape_full, x0=[d], args=(B, s1, s2, l, Vd), method='L-BFGS-B', jac=None,
        #                options={'gtol': 1e-5, 'disp': True})
        # d = res.x[0]
        print(f'Ba_res={B} d_res/diam={d/diam}')
        # X, Sigma = shape_full([d], B, s1, s2, l)
        X_psi, Sigma_psi, X_x, Sigma_x = shape_full([d], B, s1, s2, l)
        width = np.abs(X_x.min())
        height = np.abs(max(Sigma_psi.max(), Sigma_x.max()))
        plt.figure(figsize=(picScale, (height/width)*picScale))
        plt.plot(X_psi/diam, Sigma_psi/diam, '.-')
        plt.plot(X_x/diam, Sigma_x/diam, '.-')
        plt.axis('equal')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.savefig('Sigma_X.png', bbox_inches='tight')
        plt.cla()
    else:
        res = minimize(compute_shape_psi, x0=[B], args=(alpha, s1, s2, l, Vd), method='L-BFGS-B', jac=None,
                       options={'gtol': 1e-6, 'disp': True})
        B = res.x[0]
        X, Sigma = shape_psi([B], alpha, s1, s2, l)

        plt.figure(figsize=(8*picScale, 1.3*picScale))
        plt.plot(X/diam, Sigma/diam, '.')
        plt.axis('equal')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.savefig('Sigma_X.png', bbox_inches='tight')
        plt.cla()

csvPattern = "slice_t=*.csv"
xmin = -4
xmax = 4
ymax = 1.5
picScale = 4
picScale1 = 16

props = {}
props['mu1'] = 0.88e-3
props['mu2'] = 0.019e-3
props['rho1'] = 997
props['rho2'] = 1.204
props['sigma'] = 72.8e-3
props['diam'] = 0.514e-3
props['grav'] = 9.8 #variable parameter
props['alpha'] = np.pi/90  #variable parameter
props['s1'] = -1
props['s2'] = 1
props['Ub'] = 4.8 #m/s

props['Vd'] = 0.2179e-9
props['l'] = 0.51*props['diam'] #np.sqrt(props['sigma']/((props['rho1'] - props['rho2'])*props['grav'])) #???? TODO: see here

props['l_range'] = np.array([0.1, 1])*props['diam']
props['B_range'] = np.array([1, 10])/props['diam']

print(f'props={props}')

pendant_drop(props, picScale1)
a = props["a"]
print("a=", a)

# props['grav'] = dpdx/(props['rho1'] - props['rho2'])
