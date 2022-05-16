import numpy as np
import pandas as pd
import glob
import os
import re
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import minimize, rosen, rosen_der, brentq
from scipy.integrate import solve_ivp
from numpy import sin, cos, pi , exp
from scipy.integrate import odeint
from scipy import integrate

from scipy.ndimage import gaussian_filter1d
from tsmoothie.smoother import *
from scipy.interpolate import UnivariateSpline
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from scipy import interpolate
from sklearn.cluster import KMeans
SMALL_SIZE = 20
MEDIUM_SIZE = 25
BIGGER_SIZE = 30

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

def calc_logl(x,mu,sd):
    """
    Helper function to calculate log-likelihood
    """
    logl = 0
    for i in x:
        logl += np.log(norm.pdf(i, mu, sd))
    return logl

def find_optimal_k(data):
    """
    Provide a numpy array, returns index to serve as cut-off
    """
    profile_logl = []
    for q in range(1,len(data)):
        n = len(data)
        s1 = data[0:q]
        s2 = data[q:]
        mu1 = s1.mean()
        mu2 = s2.mean()
        sd1 = s1.std()
        sd2 = s2.std()
        sd_pooled = np.sqrt((((q-1)*(sd1**2)+(n-q-1)*(sd2**2)) / (n-2)))
        profile_logl.append(calc_logl(s1,mu1,sd_pooled) + calc_logl(s2,mu2,sd_pooled))
    return np.argmax(profile_logl)

def get_smooth(data, smoothing=10, return_df=False):
    if return_df:
        return pd.DataFrame(data)

    df = pd.DataFrame(data).sort_values(by=0).reset_index(drop=True).rolling(smoothing).mean().dropna()

    # first derivatives
    df['dx'] = np.gradient(df[0])
    df['dy'] = np.gradient(df[1])

    df['dx'] = df.dx.rolling(smoothing, center=True).mean()
    df['dy'] = df.dy.rolling(smoothing, center=True).mean()

    # second derivatives
    df['d2x'] = np.gradient(df.dx)
    df['d2y'] = np.gradient(df.dy)

    df['d2x'] = df.d2x.rolling(smoothing, center=True).mean()
    df['d2y'] = df.d2y.rolling(smoothing, center=True).mean()


    # calculation of curvature from the typical formula
    df['curvature'] = df.eval('abs(dx * d2y - d2x * dy) / (dx * dx + dy * dy) ** 1.5')
    # mask = curvature < 100

    df['curvature'] = df.curvature.rolling(smoothing, center=True).mean()

    df.dropna(inplace=True)
    return df

def curvature(x, y_in, a):
    x = x[::-1]
    y_in = y_in[::-1]
    #first and second derivative
    sigma = 10
    spl = UnivariateSpline(x, y_in)
    x = np.linspace(min(x), max(x), 100)
    y = spl(x)
    # smoother = ConvolutionSmoother(window_len=50, window_type='ones')
    # smoother.smooth(y_in)
    # y = smoother.smooth_data[0]
    dx_t = np.gradient(x)
    d2x_dt2 = np.gradient(dx_t)
    dy_t = np.gradient(y)
    d2y_dt2 = np.gradient(dy_t)
    curv = np.abs(d2x_dt2*dy_t - dx_t*d2y_dt2)/((dx_t)**2 + (dy_t)**2)**(3/2)
    # x1 = gaussian_filter1d(x, sigma=sigma, order=1, mode='reflect')
    # x2 = gaussian_filter1d(x1, sigma=sigma, order=1, mode='reflect')
    # y1 = gaussian_filter1d(y, sigma=sigma, order=1, mode='reflect')
    # y2 = gaussian_filter1d(y1, sigma=sigma, order=1, mode='reflect')
    # curv = np.abs(x1*y2 - y1*x2) / np.power(x1**2 + y1**2, 3./2)
    return curv

def kMeansRes(scaled_data, k, alpha_k=0.12):
    '''
    Parameters
    ----------
    scaled_data: matrix
        scaled data. rows are samples and columns are features for clustering
    k: int
        current k for applying KMeans
    alpha_k: float
        manually tuned factor that gives penalty to the number of clusters
    Returns
    -------
    scaled_inertia: float
        scaled inertia value for current k
    '''

    inertia_o = np.square((scaled_data - scaled_data.mean(axis=0))).sum()
    # fit k-means
    kmeans = KMeans(n_clusters=k, random_state=0).fit(scaled_data)
    scaled_inertia = kmeans.inertia_ / inertia_o + alpha_k * k
    return scaled_inertia


def chooseBestKforKMeans(scaled_data, k_range):
    ans = []
    for k in k_range:
        scaled_inertia = kMeansRes(scaled_data, k)
        ans.append((k, scaled_inertia))
    results = pd.DataFrame(ans, columns = ['k','Scaled Inertia']).set_index('k')
    best_k = results.idxmin()[0]
    return best_k, results

def sort_names(image_files):
    file_names = [os.path.basename(string) for string in image_files]
    times = [(re.findall("\d+\.\d+", string)[0], string) for string in file_names]
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
    # Sigma_psi = soln_psi.y[0]
    # X_psi = soln_psi.y[1]
    X_psi = soln_psi.y[0]
    Sigma_psi = soln_psi.y[1]
    return X_psi, Sigma_psi

def compute_shape_psi(B_, alpha, s1, s2, l, Vd):
    X_psi, Sigma_psi = shape_psi(B_, alpha, s1, s2, l)
    return target_fun(X_psi, Sigma_psi, Vd)

#
#w'=f''=(1 + w**2)*(1/f + np.sqrt(1 + w**2)*(s1*s2*x/l**2 - B))
#f'=w
def shape_full(d_, B, a, s1, s2, l):
    X_psi, Sigma_psi = shape_psi([B], pi/6, s1, s2, l)
    d = d_[0]
    X0 = -X_psi[-1]
    W0 = (Sigma_psi[-1] - Sigma_psi[-2])/(X_psi[-1] - X_psi[-2])
    f0 = Sigma_psi[-1]
    print(f'X0/a={X0/a} f0/a={f0/a} W0={W0} d/a={d/a}, f()={fun_x(X0, (W0, f0), s1, s2, l, B)}')
    soln_x = solve_ivp(fun_x, (X0, d), (W0, f0), method='BDF', args=(s1, s2, l, B),
                       min_step=1e-6, max_step=1e-3, rtol=1e-15, dense_output=True)
    X_x = -soln_x.t
    W_x = soln_x.y[0]
    Sigma_x = soln_x.y[1]
    X = np.concatenate([X_psi, X_x])
    Sigma = np.concatenate([Sigma_psi, Sigma_x])
    # print('X_x/a', X_x/a)
    # print('W_x', W_x)
    # print('Sigma_x/a', Sigma_x/a)
    return X_psi, Sigma_psi, X_x, Sigma_x
    # return X, Sigma

def compute_shape_full(d_, B, a, s1, s2, l, Vd):
    X_psi, Sigma_psi, X_x, Sigma_x = shape_full(d_, B, a, s1, s2, l)
    X = np.concatenate([X_psi, X_x])
    Sigma = np.concatenate([Sigma_psi, Sigma_x])
    return target_fun(X, Sigma, Vd)


# s1 = -1 for pendant drop
# s2 = 1 for rhod - rhoa > 0
def pendant_drop(props, picScale, s1 = -1, s2 = 1, N = 1000):
    Vd = props['Vd']
    alpha = props['alpha']
    a = (3*Vd/(4*np.pi))**(1/3)
    props["a"] = a
    d = a*props['d/a']
    l = np.sqrt(props['sigma']/((props['rho1'] - props['rho2'])*props['grav']))
    iBo = (l/a)**2
    B = 2/(a*(4/(2 + cos(alpha)**3 - 3*cos(alpha)))**(1/3))
    # B *= 0.9
    print(f'Ba_guess={B*a}, d/a={d/a}, a={a}, l={l/a}, a/Bo={iBo*a}')
    if True:
        # res = minimize(compute_shape_full, x0=[d], args=(B, a, s1, s2, l, Vd), method='L-BFGS-B', jac=None,
        #                options={'gtol': 1e-5, 'disp': True})
        # d = res.x[0]
        print(f'Ba_res={B*a} d_res/a={d/a}')
        # X, Sigma = shape_full([d], B, a, s1, s2, l)
        X_psi, Sigma_psi, X_x, Sigma_x = shape_full([d], B, a, s1, s2, l)
        width = np.abs(X_x.min())
        height = np.abs(max(Sigma_psi.max(), Sigma_x.max()))
        plt.figure(figsize=(picScale, (height/width)*picScale))
        plt.plot(X_psi/a, Sigma_psi/a, '.-')
        plt.plot(X_x/a, Sigma_x/a, '.-')
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
        plt.plot(X/a, Sigma/a, '.')
        plt.axis('equal')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.savefig('Sigma_X.png', bbox_inches='tight')
        plt.cla()

csvPattern = "slice_t=*.csv"
xmax = 8
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
props['grav'] = 9.8
props['alpha'] = 110*np.pi/180
props['d/a'] = 1.295828280810274

df = pd.read_csv('/Users/weugene/Desktop/toDelete/points_drop.csv', sep=',', usecols=['Points_0','Points_1'])
df = df[df['Points_1'] > 0]
left_x = df['Points_0'].min()
df['Points_0'] -= left_x
df = df.sort_values('Points_0')
df = df.reset_index(drop=True)
# print(df)
props['Vd'] = volume(df['Points_0'].values, df['Points_1'].values)*props['diam']**3#0.2179e-9
print(f'props={props}')

pendant_drop(props, picScale1, s1 = -1, s2 = 1, N = 1000)
a = props["a"]
print("a=", a)


csvnames = glob.glob(f'./{csvPattern}', recursive = False)
csvnames = sort_names(csvnames)
print('Found pvd files in:',csvnames)

for file in csvnames:
    res = pd.read_csv(file, sep = ',', usecols=['Points_0','Points_1','u.x_0','u.x_1','u.x_2','u.x_Magnitude'])
    left_x = res['Points_0'].min()
    right_x = res['Points_0'].max()
    length = right_x - left_x
    res['Points_0'] -= left_x
    plt.figure(figsize=(8*picScale, 1.3*picScale))
    x = res['Points_0'].values
    y = res['Points_1'].values
    points = np.c_[x, y]

    ind = np.abs(y) < 0.05
    # ind = (x > 1) & (np.abs(y) < 0.25)

    # Sum_of_squared_distances = np.zeros((9,2))
    K = range(1,10)
    # for i,k in enumerate(K):
    #     km = KMeans(n_clusters=k)
    #     km = km.fit(points[ind])
    #     Sum_of_squared_distances[i]=[k, km.inertia_]
    # plt.plot(Sum_of_squared_distances[:,0], Sum_of_squared_distances[:,1])

    # best_k, results = chooseBestKforKMeans(points, K)
    # plt.plot(results)
    best_k = find_optimal_k(points)
    # plt.savefig(file[:-4]+'_xz.png', bbox_inches='tight')
    # plt.cla()

    kmeans = KMeans(n_clusters=best_k, random_state=0).fit(points[ind])
    centers = sorted(kmeans.cluster_centers_, key=lambda x: x[0])
    print(f"best_k={best_k} kmeans.cluster_centers_={centers}")
    center1 = centers[-2]
    ind = np.sqrt((points[:,0] - center1[0])**2 +(points[:,1] - center1[1])**2) < 0.25

    df = get_smooth(points[ind])
    # clf = NearestNeighbors(2).fit(points)
    # G = clf.kneighbors_graph()
    # T = nx.from_scipy_sparse_matrix(G)
    #
    #
    # paths = [list(nx.dfs_preorder_nodes(T, i)) for i in range(len(points))]
    # mindist = np.inf
    # minidx = 0
    #
    # for i in range(len(points)):
    #     p = paths[i]           # order of nodes
    #     ordered = points[p]    # ordered nodes
    #     # find cost of that order by the sum of euclidean distances between points (i) and (i+1)
    #     cost = (((ordered[:-1] - ordered[1:])**2).sum(1)).sum()
    #     if cost < mindist:
    #         mindist = cost
    #         minidx = i
    # opt_order = paths[minidx]
    # x = x[opt_order]
    # y = y[opt_order]
    # print("x:", x, "y:", y)
    # dx_t = np.gradient(x)
    # d2x_dt2 = np.gradient(dx_t)
    # dy_t = np.gradient(y)
    # d2y_dt2 = np.gradient(dy_t)
    # curvature = np.abs(d2x_dt2*dy_t - dx_t*d2y_dt2)/((dx_t)**2 + (dy_t)**2)**(3/2)
    # plt.plot(points[:,0], points[:,1], 'r.')

    plt.plot(x, y, '.', ms=2)
    plt.plot(df[0], df[1], '.', ms=2)
    plt.plot([0,xmax],[-0.5,-0.5], c='0.55')
    plt.plot([0,xmax],[ 0.5, 0.5], c='0.55')
    plt.xlim(0, xmax)
    plt.ylim(-0.5, 0.5)
    plt.axis('equal')
    plt.grid(True)
    # plt.savefig(file[:-3]+'eps', bbox_inches='tight')
    plt.savefig(file[:-3]+'png', bbox_inches='tight')
    plt.cla()

    # curv = curvature(x, y, a)
    # print('curvature', curv)
    # plt.plot(curv, c='0.55')
    # plt.axis('auto')
    # # plt.show()
    # plt.savefig(file[:-4]+'_curvature.png', bbox_inches='tight')
    # plt.cla()