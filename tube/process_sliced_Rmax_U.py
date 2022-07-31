from process_sliced_bubble import *
import json
from scipy.signal import find_peaks

#find point where it start to decrease
def find_extremum(coords, find_decrease=True, eps=1e-3):
    my_list_x = coords[:,0]
    my_list_y = coords[:,1]
    my_result = len(my_list_y) - 1
    if find_decrease:
        peaks, _ = find_peaks(my_list_y, height=0, distance=100)
        print(f"find_decrease peaks={peaks}  my_list_y={my_list_y[peaks]}")
        my_result = peaks[0]
        # for index in range(1, len(my_list_y) - 1):
        #     # if my_list_x[index] < -0.01:
        #     #     continue
        #     if my_list_y[index] - my_list_y[index + 1] > eps:
        #         my_result = index
        #         break
    else:
        peaks, _ = find_peaks(-my_list_y, height=0, distance=100)
        print(f"find_increasing peaks={peaks}  my_list_y={my_list_y[peaks]}")
        my_result = peaks[0]
        # for index in range(1, len(my_list_y) - 1):
        #     # if my_list_x[index] < -0.01:
        #     #     continue
        #     if my_list_y[index + 1] - my_list_y[index] > eps:
        #         my_result = index
        #         break
    return my_result

def give_coord(x, y, df_in, index, side='up'):
    df = df_in.copy(deep=True)
    # find the start point
    # df_cluster = df[(df['label'] == index)] # take points of cluster "index"
    if side == 'up':
        df_ind = df[(df['label'] == index) & (df['y'] > 0)] # take points of cluster "index" and above y>0
    else:
        df_ind = df[(df['label'] == index) & (df['y'] <= 0)] # take points of cluster "index" and above y>0
    start = list(sorted(zip(df_ind['x'].values, df_ind['y'].values), key = lambda x: x[0])[-1]) # sort by x and take the last right element
    print(f'start:{start}')

    # x_tip = start[0]
    # x -= x_tip
    # df['x'] -= x_tip
    # df_cluster['x'] -= x_tip
    # start[0] -= x_tip

    coords = np.c_[x,y]
    if side == 'up':
        coords = coords[coords[:,1] > 0]
    else:
        coords = coords[coords[:,1] <= 0]
    clustering = DBSCAN(eps=0.02, min_samples=2).fit(coords)
    print(f'Found N={clustering} clusters')
    df_compare = pd.DataFrame({'x': coords[:,0], 'y': coords[:,1], 'label': clustering.labels_} )
    # print(f'df: {df_compare}')
    label = df_compare[(df_compare['x'] == start[0]) & (df_compare['y'] == start[1])]['label'].values[0]
    print(f'label: {label}')
    ind = df_compare['label'] == label
    coords = np.c_[df_compare['x'][ind].values, df_compare['y'][ind].values]
    coords = optimized_path(coords, start)
    return coords

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
props['Vd'] = 0.2179e-9 #estiamte volume
props["a"] = 0.000373 #estiamte radius of drop
print(f'props={props}')

csvnames = glob.glob(f'./{csvPattern}', recursive = False)
csvnames = sort_names(csvnames)
print('Found pvd files in:',csvnames)

outputs = {'t':[], 'curvature_tip': [], 'U_tip':[], 'x_tip':[], 'rmax':[]}

for ifile, file in enumerate(csvnames):
    time = get_time(file)
    # Debugging
    # if time != 0.556045:
    #     continue
    # if time > 7.61:
    #     continue
    print(f'file: {file} time: {time}')
    res = pd.read_csv(file, sep = ',', usecols=['Points_0','Points_1','u.x_0','u.x_1','u.x_2','u.x_Magnitude'])
    left_x = res['Points_0'].min()
    right_x = res['Points_0'].max()
    length = right_x - left_x
    res['Points_0'] -= left_x
    plt.figure(figsize=(8*picScale, 1.3*picScale))
    x = res['Points_0'].values
    y = res['Points_1'].values
    U = res['u.x_Magnitude'].values
    try:
        df, centers = find_df_centers(res, width=0.125) # some points |y| < width
        # draw all centers
        # for index, row in centers.iterrows():
        #     x_circle, y_circle, curvature_tip = compute_curvature(index, row, df, props["a"]))
        #     plt.plot(x_circle, y_circle, 'c.', ms=2)
        # draw the second circle from right
        index, row = list(centers.index)[-2], centers.iloc[-2]
        print('row', row, 'index', index)
        x_circle, y_circle, curvature_tip = compute_curvature(index, row, df, props["a"])

        # find the start point
        df_cluster = df[(df['label'] == index)].copy(deep=True) # take points of cluster "index"
        df_ind = df[(df['label'] == index) & (df['y'] > 0)] # take points of cluster "index" and above y>0
        # print(f'df_cluster {df_cluster}')
        start = list(sorted(zip(df_ind['x'].values, df_ind['y'].values), key = lambda x: x[0])[-1]) # sort by x and take the last right element
        # shift bubble to the beginning
        coords_up = give_coord(x, y, df, index, side='up')
        coords_down = give_coord(x, y, df, index, side='down')
    except:
        try:
            df, centers = find_df_centers(res, width=0.1) # some points |y| < width
            # draw all centers
            # for index, row in centers.iterrows():
            #     x_circle, y_circle, curvature_tip = compute_curvature(index, row, df, props["a"]))
            #     plt.plot(x_circle, y_circle, 'c.', ms=2)
            # draw the second circle from right
            index, row = list(centers.index)[-2], centers.iloc[-2]
            print('row', row, 'index', index)
            x_circle, y_circle, curvature_tip = compute_curvature(index, row, df, props["a"])

            # find the start point
            df_cluster = df[(df['label'] == index)].copy(deep=True) # take points of cluster "index"
            df_ind = df[(df['label'] == index) & (df['y'] > 0)] # take points of cluster "index" and above y>0
            # print(f'df_cluster {df_cluster}')
            start = list(sorted(zip(df_ind['x'].values, df_ind['y'].values), key = lambda x: x[0])[-1]) # sort by x and take the last right element
            # shift bubble to the beginning
            coords_up = give_coord(x, y, df, index, side='up')
            coords_down = give_coord(x, y, df, index, side='down')
        except:
            print('Can\'t find start point')
            continue
    print(f'start:{start}')

    x_tip = start[0]
    U_tip = df[(df['x'] == start[0]) & (df['y'] == start[1])]['U'].values[0]
    x_circle -= x_tip
    x -= x_tip
    df['x'] -= x_tip
    df_cluster['x'] -= x_tip
    start[0] -= x_tip
    coords_up[:,0] -= x_tip
    coords_down[:,0] -= x_tip
    i_rmax_up = min(find_extremum(coords_up, True), np.argmin(coords_up[:,0])) # choose index between point where y begins decreasing and the tail of the bubble
    rmax_up = coords_up[i_rmax_up,1]
    i_rmax_down = min(find_extremum(coords_down, False), np.argmin(coords_down[:,0])) # choose index between point where y begins decreasing and the tail of the bubble
    rmax_down = coords_down[i_rmax_down,1]
    print(f"up: {find_extremum(coords_up, True)} {np.argmin(coords_up[:,0])}")
    print(f"down: {find_extremum(coords_down, False)} {np.argmin(coords_down[:,0])}")
    print(f'rmax_up[{i_rmax_up}]={rmax_up} rmax_down[{i_rmax_down}]={rmax_down}')
    # Output to the file
    outputs['t'].append(time)
    outputs['curvature_tip'].append(curvature_tip)
    outputs['U_tip'].append(U_tip)
    outputs['x_tip'].append(x_tip)
    outputs['rmax'].append(max(np.abs(rmax_down), np.abs(rmax_up)))

    # Draw
    df = df.values
    df_cluster = df_cluster.values

    width = np.abs(xmax - xmin)
    height = 1
    plt.figure(figsize=(picScale1, (height/width)*picScale1))
    plt.plot(x_circle, y_circle, 'c-', ms=2)
    # plt.plot(coords_up[::5,0], coords_up[::5,1], '-', lw=4)
    # plt.plot(coords_down[::5,0], coords_down[::5,1], '-', lw=4)
    plt.plot(x, y, '.', ms=2) #all points
    plt.plot(df_cluster[:,0], df_cluster[:,1], 'r.', ms=2) #chosen only 1 cluster
    # plt.plot(df[:,0], df[:,1], 'r.', ms=2) #chosen for clustering |y| <0.1
    # if rmax_up > rmax_down:
    plt.plot([coords_up[i_rmax_up, 0]], [coords_up[i_rmax_up, 1]], '.', ms=5) #Rmax point
    # else:
    plt.plot([coords_down[i_rmax_down, 0]], [coords_down[i_rmax_down, 1]], '.', ms=5) #Rmax point
    plt.plot([xmin,xmax],[-0.5,-0.5], c='0.55')
    plt.plot([xmin,xmax],[ 0.5, 0.5], c='0.55')
    plt.xlim(xmin, xmax)
    plt.ylim(-0.5, 0.5)
    # displaying the title
    plt.title(f"$t={time}$")
    # plt.axis('equal')
    plt.grid(True)
    plt.savefig(file[:-3]+'eps', bbox_inches='tight')
    plt.savefig(file[:-3]+'png', bbox_inches='tight')
    plt.cla()

with open("output_sliced_Rmax_U.json", "w") as f:
    json.dump(outputs, f)