import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
import glob, os, re, json

def sort_names(image_files):
    file_names = [os.path.basename(string) for string in image_files]
    times = [(float(re.findall("\d+\.\d+", string)[0]), string) for string in file_names]
    times = sorted(times, key=lambda x: x[0])
    print(f"Time: {times[0][0]} -- {times[-1][0]}")
    image_files = [t[1] for t in times]
    return image_files

def give_time(file):
    filename = os.path.basename(file)
    time = re.findall("\d+\.\d+", filename)[0]
    return time

def plot_graph(list_x, list_y, names, xtitle, ytitle, image_name, list_x_fill=[], list_y_fill=[], mode=[], \
               dash=['solid', 'dot', 'dash', 'longdash'], \
               colors=['blue', 'red', 'hsv(120,100,100)', 'green', 'black' ], \
               marker_size=15, xrange =[], yrange = [], \
               marker_style = ['circle', 'triangle-up', 'triangle-down','square', 'diamond', 'cross',  'x-thin', 'cross-thin' ], \
               width=1000, height=500, path='./', yanchor='center', y0_anchor=0.01, xanchor='left', x0_anchor=0.3, logx=False, logy=False):
    if mode == []:
        for i in range(len(list_x)):
            mode.append('lines+markers')


    while len(marker_style) < len(list_x):
        marker_style[:] = marker_style[:] + marker_style[:]
    figborderlinesize = 0.7
    legborderlinesize = 0.7
    yaxis = dict(
        tickfont = dict(
            family = 'Times New Roman',
            size = 20,
            color = 'black'
        ),
        titlefont = dict(
            family = 'Times New Roman',
            size = 25,
            color = 'black'
        )
    )
    if logy:
        yaxis['type'] = "log"
    xaxis = dict(
        tickfont = dict(
            family = 'Times New Roman',
            size = 20,
            color = 'black'
        ),
        titlefont = dict(
            family = 'Times New Roman',
            size = 25,
            color = 'black'
        )
    )
    if logx:
        xaxis['type'] = "log"

    axis_style = dict(showline=True, gridwidth=1, gridcolor='lightgrey', linewidth=figborderlinesize, linecolor='black', mirror=True, ticks='outside', tickfont = dict(family = 'Times New Roman', size = 20, color = 'black'))
    bg_style = {'plot_bgcolor': 'rgba(255, 255, 255, 1)', 'paper_bgcolor': 'rgba(255, 255, 255, 1)',}


    fig = go.Figure()
    k = len(list_x)
    n_fill = len(list_x_fill)
    if len(list_x_fill) == 2 and len(list_y_fill) == 2:
        fig.add_trace(go.Scatter(x=list_x_fill[1], y=list_y_fill[1], name=names[k+1], mode='lines', fillcolor='blueviolet', line_color='blueviolet', fill='tozeroy')) # fill to trace0 y
        fig.add_trace(go.Scatter(x=list_x_fill[0], y=list_y_fill[0], name=names[k], mode='lines', fillcolor='lightsteelblue',     line_color='indigo', fill='tozeroy')) # fill down to xaxis
    for i,x in enumerate(list_x):
        print('Plot curve number:', i)
        y = np.asarray(list_y[i])
        fig.add_trace(go.Scatter(x=x, y=y, name=names[i],
                                 mode=mode[i],
                                 marker=dict(
                                     size=marker_size,
                                     line=dict(width=1)
                                 ),
                                 marker_symbol=marker_style[i],
                                 line=dict(width=2, dash=dash[i]),
                                 textfont=dict(
                                     family="Times New Roman",
                                     size=18,
                                     color="LightSeaGreen")
                                 ))
        if colors != []:
            fig['data'][i + n_fill]['marker']['line']['color'] = colors[i]
            fig['data'][i + n_fill]['line']['color'] = colors[i]
    fig.update_xaxes(axis_style)
    fig.update_yaxes(axis_style)

    fig.update_layout(
        width = width,
        height = height,
        xaxis_title=xtitle,
        yaxis_title=ytitle,
        yaxis = yaxis,
        xaxis = xaxis,
        showlegend=True
    )
    fig.update_layout(bg_style)

    fig.update_layout(legend=dict(
        bgcolor="White",
        bordercolor="Black",
        borderwidth=figborderlinesize
    ))
    fig.update_layout(font=dict(
        family="Times New Roman",
        size=20,
        color="Black"
    ))
    fig.update_layout(
        autosize=False,
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0.1
        ),
        #     paper_bgcolor="LightSteelBlue",
    )
    fig.update_layout(legend=dict(
        yanchor=yanchor,
        y=y0_anchor,
        xanchor=xanchor,
        x=x0_anchor
    ))
    if len(xrange) == 2:
        fig.update_xaxes(range=xrange)
    if len(yrange) == 2:
        fig.update_yaxes(range=yrange)

    # fig.update_xaxes(type="log")
    # fig.update_yaxes(type="log")

    #fig.show()
    fn = path + image_name
    print('Write image to file:', fn)
    fig.write_image(str(Path(fn)))
    # fig.write_image(str(Path(fn)), engine="kaleido")
    print("Successfully generated:", fn)

if __name__ == "__main__":
    path = '/Users/weugene/basilisk/work/tube/res27/'
    with open(path + "output_sliced_Rmax_U.json", "r") as f:
        outputs = json.load(f)
    df = pd.read_csv(path + 'for_excel_table.txt', sep=' |\t', header=0)
    tmin = df['t'].values.min()
    tmax = df['t'].values.max()
    fn = "x_nose_mean_tip.png"
    plot_graph([df['t'].values, df['t'].values, outputs['t']], [df['x_nose'].values, df['x_mean'].values, outputs['x_tip']], \
               ['X nose', 'X mean', 'X tip'], \
               dash=['solid', 'dot', 'dash', 'dot'], \
               xtitle="t", ytitle="x", image_name=fn[:-3]+'pdf', mode=['lines', 'lines', 'lines', 'lines'], \
               colors=['red', 'blue', 'black' ], \
               marker_size=1, width=1000, height=500, path=path, yanchor='bottom', y0_anchor=0.01, xanchor='right', x0_anchor=0.99)

    fn = "u_nose_mean_estim_tip.png"
    umean_estim = np.diff(df['x_mean'].values)/np.diff(df['t'].values)
    plot_graph([df['t'].values, df['t'].values, outputs['t']], [np.diff(df['x_nose'].values)/np.diff(df['t'].values), umean_estim, outputs['U_tip']], \
               ['U nose', 'U mean', 'U tip'], \
               dash=['solid', 'dot', 'dash'], \
               xtitle="t", ytitle="U", image_name=fn[:-3]+'pdf', mode=['lines', 'lines', 'lines'], \
               colors=['red', 'blue', 'black' ], \
               marker_size=1, width=1000, height=500, path=path, yanchor='bottom', y0_anchor=0.01, xanchor='left', x0_anchor=0.3)


    csvPattern = "slice_t=*.csv"
    csvnames = glob.glob(f'{path}/{csvPattern}', recursive = False)
    csvnames = sort_names(csvnames)
    print('Found pvd files in:',csvnames)

    list_t = []
    list_unose = []
    # list_umean = []

    for ifile, file in enumerate(csvnames):
        time = give_time(file)
        print(f'file: {file}, time={time}')
        res = pd.read_csv(path + file, sep = ',', usecols=['Points_0','Points_1','u.x_0','u.x_1','u.x_2','u.x_Magnitude'])
        row = res.iloc[res['Points_0'].idxmax()]
        list_unose.append(row['u.x_0'])
        list_t.append(float(time))

    fn = "u_nose_mean_tip_markers.png"
    plot_graph([list_t, df['t'].values, outputs['t']], [list_unose, df['UmeanV'].values, outputs['U_tip']], \
               ['U nose', 'U mean', 'U tip'], \
               dash=['solid', 'dot', 'dash'], \
               xtitle="t", ytitle="U", image_name=fn[:-3]+'pdf', mode=['lines+markers', 'lines+markers', 'lines+markers'], \
               colors=['red', 'black', 'black' ], \
               marker_size=5, width=1000, height=500, path=path, yanchor='top', y0_anchor=0.99, xanchor='right', x0_anchor=0.99, logx=False, logy=False)

    # list_unose=np.asarray(list_unose)
    # fn = "u_nose_markers.png"
    # plot_graph([list_t, outputs['t']], [list_unose - 1.6, outputs['U_tip']], \
    #            ['$U_{\\text{nose}} - \\hat{U}$'], \
    #            dash=['solid', 'dot', 'dot'], \
    #            xtitle="t", ytitle='$U_{\\text{nose}} - \\hat{U}$', image_name=fn[:-3]+'pdf', mode=['lines+markers', 'lines+markers', 'lines+markers'], \
    #            colors=['red', 'black', 'black' ], \
    #            marker_size=5, width=1000, height=500, path=path, yanchor='top', y0_anchor=0.99, xanchor='right', x0_anchor=0.99, logx=False, logy=False)


    # print(list_t, list_unose)
    fn = "u_nose_mean_markers_log.png"
    plot_graph([list_t, df['t'].values, outputs['t']], [list_unose, df['UmeanV'].values, outputs['U_tip']], \
               ['U nose', 'U mean', 'U tip'], \
               dash=['solid', 'dot', 'dash'], \
               xtitle="t", ytitle="U", image_name=fn[:-3]+'pdf', mode=['lines+markers', 'lines+markers', 'lines+markers'], \
               colors=['red', 'black', 'black' ], \
               marker_size=5, width=1000, height=500, path=path, yanchor='top', y0_anchor=0.99, xanchor='right', x0_anchor=0.99, logx=True, logy=True)

    # list_unose=np.asarray(list_unose)
    # fn = "u_nose_markers_log.png"
    # plot_graph([list_t], [list_unose - 1.6], \
    #            ['$U_{\\text{nose}} - \\hat{U}$'], \
    #            dash=['solid', 'dot', 'dot'], \
    #            xtitle="t", ytitle='$U_{\\text{nose}} - \\hat{U}$', image_name=fn[:-3]+'pdf', mode=['lines+markers', 'lines+markers', 'lines+markers'], \
    #            colors=['red', 'black', 'black' ], \
    #            marker_size=5, width=1000, height=500, path=path, yanchor='top', y0_anchor=0.99, xanchor='right', x0_anchor=0.99, logx=True, logy=True)
