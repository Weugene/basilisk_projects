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

def func_curvature(pars, x, y):
    xc, yc, R = pars
    return sum([np.abs((xi - xc)**2 + (yi - yc)**2 - R**2) for xi, yi in zip(x, y)])

def find_df_centers(res):
    x = res['Points_0'].values
    y = res['Points_1'].values
    points = np.c_[x, y]

    ind = np.abs(y) < 0.1
    clustering = DBSCAN(eps=0.01, min_samples=2).fit(points[ind])
    df = pd.DataFrame({'x': points[ind][:,0], 'y': points[ind][:,1], 'label': clustering.labels_} )
    centers = df.groupby('label').mean() #label becomes as an index
    centers.sort_values(by=['x'], inplace=True)
    print(f'centers: {centers}')
    return df, centers

def compute_curvature(index, row, df, a):
    center_x, center_y = row['x'], row['y']
    df_label = df[df['label'] == index]
    res = minimize(func_curvature, x0=[center_x, 0, 25*a], args=(df_label['x'].values, df_label['y'].values), method='L-BFGS-B', jac=None,
                   options={'gtol': 1e-6, 'disp': True})
    print(res.x)
    curvature = 2/res.x[2]
    print(f'RES: center_x={res.x[0]}, center_y={res.x[1]}, R={res.x[2]}, curvature={curvature}')
    x_circle, y_circle = res.x[0] + res.x[2]*np.cos(phi), res.x[1] + res.x[2]*np.sin(phi)
    return x_circle, y_circle, curvature

def distance(P1, P2):
    """
    This function computes the distance between 2 points defined by
     P1 = (x1,y1) and P2 = (x2,y2)
    """
    return ((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2) ** 0.5

def optimized_path(coords, start=None):
    """
    This function finds the nearest point to a point
    coords should be a list in this format coords = [ [x1, y1], [x2, y2] , ...]

    """
    if type(coords) != list:
        coords = [[x[0], x[1]] for x in coords]

    if start is None:
        start = coords[0]
    pass_by = coords
    path = [start]
    pass_by.remove(start)
    while pass_by:
        nearest = min(pass_by, key=lambda x: distance(path[-1], x))
        path.append(nearest)
        pass_by.remove(nearest)
    return np.asarray(path)

def uniform_interpolation(x, y, xn, yn):
    N = 1000
    print(x.min(), xn.min())
    xmin = max(x.min(), xn.min())
    xmax = min(x.max(), xn.max())
    xx = np.linspace(xmin, xmax, N)
    f = interp1d(x, y)
    yy = f(xx)

    f = interp1d(xn, yn)
    yyn = f(xx)
    # return LA.norm(yy - yyn, ord=1)*(xmax - xmin)
    return np.abs(yy - yyn).sum()

# `values` should be sorted
def get_closest_ind(array, values):
    # make sure array is a numpy array
    array = np.array(array)

    # get insert positions
    idxs = np.searchsorted(array, values, side="left")

    # find indexes where previous index is closer
    prev_idx_is_less = ((idxs == len(array))|(np.fabs(values - array[np.maximum(idxs-1, 0)]) < np.fabs(values - array[np.minimum(idxs, len(array)-1)])))


    try:
        idxs[prev_idx_is_less] -= 1
        return idxs[0]
    except:
        return idxs

def fit_curve_err(par, props, coords):
    tail_y, B0, l = par
    diam = props['diam']
    d = diam*props['d/diam'] # TODO: be careful
    global iter
    print(f'fit_curve_err: par: {par}, d/diam: {d/diam}')

    coords_filtered = coords[(coords[:,0] > -np.abs(d)/diam) & (coords[:,0] < 0)] # take points until d

    id = get_closest_ind(coords_filtered[:,0], -np.abs(d)/diam)
    point_left = coords_filtered[id]
    print(f'Closest left point/diam: {point_left/diam}')
    # print(f'coords: {coords}')
    # print(f'coords_filtered: {coords_filtered}')

    # B = props['B']
    s1 = props['s1']
    s2 = props['s2']
    X_psi, Sigma_psi, X_x, Sigma_x = shape_full([np.abs(d)], B0, s1, s2, l)
    # width = np.abs(X_x.min())
    # height = np.abs(max(Sigma_psi.max(), Sigma_x.max()))
    # plt.figure(figsize=(picScale, (height/width)*picScale))
    # plt.plot(X_psi/diam, Sigma_psi/diam)
    # plt.plot(X_x/diam, Sigma_x/diam)
    # plt.grid(True)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.axis('equal')
    # plt.savefig(f'fit_curve_err_{iter}.png', bbox_inches='tight')
    # plt.cla()

    iter += 1
    # print(X_psi, Sigma_psi, X_x, Sigma_x)
    print(f'X_psi: {X_psi.min()} {X_psi.max()}')
    print(f'X_x: {X_x.min()} {X_x.max()}')
    coords_filtered1 = coords_filtered[coords_filtered[:,0] > X_psi.min()/diam]
    coords_filtered2 = coords_filtered[coords_filtered[:,0] <= X_psi.min()/diam]
    # print(f'coords_filtered1: {coords_filtered1}')
    # print(f'coords_filtered2: {coords_filtered2}')
    err1 = err2 = 0
    if len(coords_filtered1):
        err1 = uniform_interpolation(X_psi/diam, Sigma_psi/diam, coords_filtered1[:,0], coords_filtered1[:,1])
    if len(coords_filtered2):
        err2 = uniform_interpolation(X_x/diam,   Sigma_x/diam,   coords_filtered2[:,0], coords_filtered2[:,1])
    return err1 + err2 + 10*distance(point_left, [props['d/diam'], tail_y/diam])

def fit_curve(par0, props, coords):
    res = dict({'x': par0})
    res = minimize(fit_curve_err, x0=par0, args=(props, coords), method='Nelder-Mead'
                   )
                   # options={'gtol': 1e-4, 'disp': False})
    d = props['d']
    print(res['x'])
    tail_y, B0, l = res['x']
    X_psi, Sigma_psi, X_x, Sigma_x = shape_full([np.abs(d)], B0, props['s1'], props['s2'], l)
    props['d'] = np.abs(d)
    props['d/diam'] = np.abs(d)/props['diam']
    props['l'] = l
    diam = props['diam']
    plt.plot(X_psi/diam, Sigma_psi/diam, '.-')
    plt.plot(X_x/diam, Sigma_x/diam, '.-')
    plt.grid(True)
    plt.savefig('shape_full.eps', bbox_inches='tight')
    plt.cla()
    return X_psi/diam, Sigma_psi/diam, X_x/diam, Sigma_x/diam, props


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
props['alpha'] = 110*np.pi/180  #variable parameter
props['d/diam'] = 1.295828280810274  #variable parameter
props['s1'] = -1
props['s2'] = 1
df = pd.read_csv('/Users/weugene/Desktop/toDelete/points_drop.csv', sep=',', usecols=['Points_0','Points_1'])
df = df[df['Points_1'] > 0]
left_x = df['Points_0'].min()
df['Points_0'] -= left_x
df = df.sort_values('Points_0')
df = df.reset_index(drop=True)
# print(df)
props['Vd'] = volume(df['Points_0'].values, df['Points_1'].values)*props['diam']**3 #0.2179e-9
print(f'props={props}')

pendant_drop(props, picScale1)
a = props["a"]
print("a=", a)


csvnames = glob.glob(f'./{csvPattern}', recursive = False)
csvnames = sort_names(csvnames)
print('Found pvd files in:',csvnames)

dpdx = (0.428276-0.510173)/(3.95884-5.19025)
# props['grav'] = dpdx/(props['rho1'] - props['rho2'])
props['l'] = np.sqrt(props['sigma']/((props['rho1'] - props['rho2'])*props['grav'])) #???? TODO: see here

phi = np.linspace(0, 2*np.pi, 100)
for ifile, file in enumerate(csvnames):
    if ifile != 2:
        continue
    print(f'file: {file}')
    res = pd.read_csv(file, sep = ',', usecols=['Points_0','Points_1','u.x_0','u.x_1','u.x_2','u.x_Magnitude'])
    left_x = res['Points_0'].min()
    right_x = res['Points_0'].max()
    length = right_x - left_x
    res['Points_0'] -= left_x
    plt.figure(figsize=(8*picScale, 1.3*picScale))
    x = res['Points_0'].values
    y = res['Points_1'].values
    df, centers = find_df_centers(res)
    # draw all centers
    # for index, row in centers.iterrows():
    #     x_circle, y_circle, curvature = compute_curvature(index, row, df, a)
    #     plt.plot(x_circle, y_circle, 'c.', ms=2)
    # draw the second circle from right
    index, row = list(centers.index)[-2], centers.iloc[-2]
    print('row', row, 'index', index)
    x_circle, y_circle, curvature = compute_curvature(index, row, df, a)

    # find the start point
    df_cluster = df[(df['label'] == index)]
    df_ind = df[(df['label'] == index) & (df['y'] > 0)]
    start = list(sorted(zip(df_ind['x'].values, df_ind['y'].values), key = lambda x: x[0])[-1])
    print(f'start:{start}')
    # shift bubble to the beginning
    x_tip = start[0]
    x_circle -= x_tip
    x -= x_tip
    df['x'] -= x_tip
    df_cluster['x'] -= x_tip
    start[0] -= x_tip

    coords = np.c_[x,y]
    coords = coords[coords[:,1] > 0]
    clustering = DBSCAN(eps=0.02, min_samples=2).fit(coords)
    print(f'Found N={clustering} clusters')
    df_compare = pd.DataFrame({'x': coords[:,0], 'y': coords[:,1], 'label': clustering.labels_} )
    # print(f'df: {df_compare}')
    label = df_compare[(df_compare['x'] == start[0]) & (df_compare['y'] == start[1])]['label'].values[0]
    print(f'label: {label}')
    ind = df_compare['label'] == label
    coords = np.c_[df_compare['x'][ind].values, df_compare['y'][ind].values]
    coords = optimized_path(coords, start)
    # print(f'coords = {coords}')
    # Find

    props['d/diam'] = 1#1.3#0.8801273617937209  #variable parameter 2.2
    props['tail_y/diam'] = 0.2
    props['l'] = 2.54107886e-04#0.00020038 #0.00020039#0.00019
    # props['Vd'] = 8e-1 #3.642e-11
    d0 = props['d/diam']*props['diam']
    tail_y = props['tail_y/diam']*props['diam']
    B0 = curvature/props['diam']
    l0 = props['l']
    props['B'] = curvature/props['diam']
    props['d'] = props['d/diam']*props['diam']
    X_psi_theor, Sigma_psi_theor, X_x_theor, Sigma_x_theor, props = fit_curve((tail_y, B0, l0), props, coords)
    # print([a for a in zip(X_psi_theor, Sigma_psi_theor)])
    # print([a for a in zip(X_x_theor, Sigma_x_theor)])
    a = props['a']
    print(f"props[d/diam] {props['d/diam']}")

    df = df.values
    df_cluster = df_cluster.values

    width = np.abs(xmax - xmin)
    height = 1
    plt.figure(figsize=(picScale1, (height/width)*picScale1))
    plt.plot(x_circle, y_circle, 'c.', ms=2)
    # plt.plot(coords[::5,0], coords[::5,1], '-', lw=4)
    plt.plot(x, y, '.', ms=2) #all points
    plt.plot(X_psi_theor, Sigma_psi_theor, 'y-')
    plt.plot(X_x_theor, Sigma_x_theor, 'y-')
    plt.plot(X_psi_theor, -Sigma_psi_theor, 'y-')
    plt.plot(X_x_theor, -Sigma_x_theor, 'y-')
    plt.plot(df_cluster[:,0], df_cluster[:,1], 'r.', ms=2) #chosen only 1 cluster
    # plt.plot(df[:,0], df[:,1], 'r.', ms=2) #chosen for clustering |y| <0.1

    plt.plot([xmin,xmax],[-0.5,-0.5], c='0.55')
    plt.plot([xmin,xmax],[ 0.5, 0.5], c='0.55')
    plt.xlim(xmin, xmax)
    plt.ylim(-0.5, 0.5)
    # plt.axis('equal')
    plt.grid(True)
    plt.savefig(file[:-3]+'eps', bbox_inches='tight')
    plt.savefig(file[:-3]+'png', bbox_inches='tight')
    plt.cla()