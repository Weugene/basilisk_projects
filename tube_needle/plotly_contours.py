########################### PLOTLY #####################################################################
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import iplot, init_notebook_mode
import plotly.io as pio
import glob, os, sys
import logging
from sys import argv

#Current PATH reading
path = os.path.abspath(os.getcwd())
print("Current PATH=" + path)

files_csv = []
files_csv += glob.glob(path + "/*.csv")
print(files_csv)
df = []


fig = go.Figure()

for i, file in enumerate(files_csv):
    df.append(pd.read_csv(file))
    mdf = df[i]
    N = mdf.shape[0]
    dfdown = mdf[(mdf['Points_1']<0)].sort_values(by=['Points_0', 'Points_1'])
    dfup = mdf[(mdf['Points_1']>0)].sort_values(by=['Points_0', 'Points_1'], ascending=False)
    mdf = pd.concat([dfdown, dfup, dfdown.head(1)])
    x = mdf['Points_0']
    y = mdf['Points_1']
    print(i)
    fig.add_trace(go.Scatter(x=x, y=y, name='t='+str(i),
                             mode='lines',
                             textfont=dict(
                                family="sans serif",
                                size=18,
                                color="LightSeaGreen")
                             ))
    fig1 = go.Figure(go.Scatter(x=x, y=y, name='t='+str(i),
                             mode='lines',
                             textfont=dict(
                                family="sans serif",
                                size=18,
                                color="LightSeaGreen")
                             ))
    fig1.update_layout(
        width = 1000,
        height = 500,
        xaxis_title='x',
        yaxis_title='y'
    )

    fig1.update_yaxes(range=[-0.5,0.5])
    iplot(fig1)
    fn = path + 'bubble_evolution_t=' + str(i) + '.pdf'
    pio.write_image(fig1, fn)
    print('File=' + fn + ' generated succesfully')

fig.update_layout(
    width = 1000,
    height = 500,
    xaxis_title='x',
    yaxis_title='y'
)



fig.update_yaxes(range=[-0.5,0.5])
iplot(fig)
fn = path + 'bubble_evolution.pdf'
pio.write_image(fig, fn)
print('File=' + fn + ' generated succesfully')
# fig.show()

