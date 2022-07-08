from process_sliced_bubble import *


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

outputs = []
phi = np.linspace(0, 2*np.pi, 100)
for ifile, file in enumerate(csvnames):
    time = get_time(file)
    # if ifile != 2:
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
    df, centers = find_df_centers(res)
    # draw all centers
    # for index, row in centers.iterrows():
    #     x_circle, y_circle, curvature = compute_curvature(index, row, df, props["a"]))
    #     plt.plot(x_circle, y_circle, 'c.', ms=2)
    # draw the second circle from right
    index, row = list(centers.index)[-2], centers.iloc[-2]
    print('row', row, 'index', index)
    x_circle, y_circle, curvature = compute_curvature(index, row, df, props["a"])

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
    # Draw

    df = df.values
    df_cluster = df_cluster.values

    width = np.abs(xmax - xmin)
    height = 1
    plt.figure(figsize=(picScale1, (height/width)*picScale1))
    plt.plot(x_circle, y_circle, 'c.', ms=2)
    # plt.plot(coords[::5,0], coords[::5,1], '-', lw=4)
    plt.plot(x, y, '.', ms=2) #all points
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