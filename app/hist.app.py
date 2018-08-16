import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from textwrap import dedent as d

import os
import json
import base64
from matplotlib.colors import TABLEAU_COLORS as tab_colors

app = dash.Dash()
asset_directory = os.getcwd() + '/assets/'
df = pd.read_csv("assets/vizapp.csv")
starting_round = 1.2
dff = df.loc[df['Round'] == starting_round]
cols = []
for item in df.columns:
    if "0" in item:
        cols.append(item)
    elif "1" in item:
        cols.append(item)
available_indicators = cols
xcolumn = available_indicators[0]
ycolumn = available_indicators[1]
def make_color_dic(dff):
    tab_keys = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    color_dic = {}
    for index, category in enumerate(dff["category"].unique()):
        color_dic[category] = tab_keys[index]
    categories = dff["category"]
    colors = []
    for i in range(len(categories)):
        colors.append(color_dic[list(categories)[i]])
    return colors

#def score_histogram(
#        dff,
#        plot_type='bar'):
#    data = [dict(
#        x = dff[dff]
#    )]
#)


def scatter_plot(
        dff,
        xcolumn,
        ycolumn,
        size=dff["Cation Heavy Atoms"],
        plot_type='scatter'):
    data = [dict(
        x=dff[dff['category'] == i][xcolumn],
        y=dff[dff['category'] == i][ycolumn],
        xlabel=xcolumn,
        ylabel=ycolumn,
        customdata=dff[dff['category'] == i]['category'],
        smiles=dff['Salt Smiles'],
        score=dff['Tanimoto Similarity Score'],
        atoms=dff['Cation Heavy Atoms'],
        relative=dff['Molecular Relative'],
        anion=dff['Anion'],
        category=dff['category'],
        name=i,
        #legendgroup=name,
        mode='markers',
        marker=dict(
            sizeref=45,
            sizemode='diameter',
            opacity=0.5,
            size=size*60,
        ),
        type=plot_type,
    ) for i in dff.category.unique()
    ]

    layout = dict(
        margin=dict(r=15, t=15),#, l=40, b=40),
        xaxis=dict(
            title=xcolumn,
        ),
        yaxis=dict(
            title=ycolumn,
        ),
        showlegend=True,
        #legend=name
    )

    return dict(data=data, layout=layout)


###for side panel

smiles = df.loc[df['Round'] == starting_round]['Salt Smiles'].iloc[1]
score = df.loc[df['Round'] == starting_round]['Tanimoto Similarity Score'].iloc[1]
atoms = df.loc[df['Round'] == starting_round]['Cation Heavy Atoms'].iloc[1]
relative =  df.loc[df['Round'] == starting_round]['Molecular Relative'].iloc[1]
anion = df.loc[df['Round'] == starting_round]['Anion'].iloc[1]
category = df.loc[df['Round'] == starting_round]['category'].iloc[1]
image_ID = df.loc[df['Round'] == starting_round]['Salt Smiles'].iloc[1].split(".")[0]


logo_filename = asset_directory + 'gains.png'
image_filename = asset_directory + image_ID + '.png'
encoded_logo = base64.b64encode(open(logo_filename, 'rb').read())
encoded_image = base64.b64encode(open(image_filename, 'rb').read())
FIGURE = scatter_plot(dff, xcolumn, ycolumn)

#styles = {
#    'pre': {
#        'border': 'thin lightgrey solid',
#        'overflowX': 'scroll'
#    }
#}
app.layout = html.Div([
    # Row 1: Header and Intro text

    html.Div([

        html.Img(src='data:image/png;base64,{}'.format(encoded_logo.decode()),
        style = {
                    'height': '100px',
                    'float': 'right',
                    'position': 'relative',
                    'bottom': '0px',
                    'left': '0px'
                },
                 ),
        html.H2('GAINS UI',
                style={
                    'position': 'relative',
                    'top': '0px',
                    'left': '10px',
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '6.0rem',
                    'color': '#4D637F'
                }),
        html.H2('for',
                style={
                    'position': 'relative',
                    'top': '0px',
                    'left': '20px',
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '4.0rem',
                    'color': '#4D637F'
                }),
        html.H2('Rapid QSPR Analysis',
                style={
                    'position': 'relative',
                    'top': '0px',
                    'left': '27px',
                    'font-family': 'Dosis',
                    'display': 'inline',
                    'font-size': '6.0rem',
                    'color': '#4D637F'
                }),

    ], className='row twelve columns',
        style={'position': 'relative', 'right': '15px'}),

    html.Div([
        html.Div([
            html.Div([
                html.P(
                    'HOVER over a salt in the graph to the left to see its structure to the right.'),
                html.P(
                    'SELECT property parameters in the dropdown menus to set graph axes.')
            ], style={'margin-left': '10px'}),
        ], className='twelve columns')

    ], className='row'),

    # Row 2: Hover Panel and Graph
    html.Div([
        html.Div([
            html.Img(id='chem_img', src='data:image/png;base64,{}'.format(encoded_image.decode()),
            style = {
                        'height': 'auto',
                        'width': '100%',
                        'float': 'top',
                        'position': 'relative',
                        'bottom': '0px',
                    }),
            html.Br(),

            html.P(
                'Description',
                style = {
                    'left': '10px'
                }
            ),
            html.P('Heavy atoms:  {}'.format(atoms),
                   id='chem_atoms',
                   style=dict(maxHeight='400px', fontSize='12px', left='10px')),
            html.P('Anion:  {}'.format(anion),
                   id='chem_anion',
                   style=dict(maxHeight='400px', fontSize='12px', left='10px')),
            html.P('Closest relative:  {}'.format(relative),
                   id='chem_relative',
                   style=dict(maxHeight='400px', fontSize='12px', left='10px')),
            html.P('Similarity score:  {:.3}'.format(score),
                   id='chem_score',
                   style=dict(maxHeight='400px', fontSize='12px', left='10px')),
            html.P('SMILES:  {}'.format(smiles.split(".")[0]),
                   id='chem_smiles',
                   style=dict(maxHeight='400px', fontSize='12px', left='10px')),
#             dcc.Markdown(d("""
#                 **Hover Data**
#
#                 Mouse over values in the graph.
#             """)),
#             html.Pre(id='hover-data')


        ], className='three columns',
            style={'margin-left': '10px'}),
        #], className='three columns', style=dict(height='300px', left='10px')),
        html.Div([
            dcc.Dropdown(
                id='xaxis-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='prediction 0'
            ),

        ], className='three columns'),

        html.Div([
            dcc.Dropdown(
                id='yaxis-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='prediction 1'
            ),

        ], className='three columns'),
    ]),
    html.Div([
        dcc.Graph(id='indicator-graphic',
        hoverData={'points': [{'customdata': 'Imidazolium'}]},
        figure=FIGURE),


        dcc.Slider(
            id='year--slider',
            min=df['Round'].min(),
            max=df['Round'].max(),
            value=df['Round'].max(),
            step=None,
            marks={str(year): str(year) for year in df['Round'].unique()}
        ),
    ], className='six columns'
    ),
    html.Div([
        dcc.Graph(id='score-histogram',
                  figure={
                      'data': [
                          go.Histogram(
                              x=dff["Cation Heavy Atoms"],
                            xbins=dict(
                               start=dff["Cation Heavy Atoms"].min(),
                               end=dff["Cation Heavy Atoms"].max(),
                               size=1,
                                ),
                          )
                      ]
                  })
    ], className='two columns'
    ),
])

# @app.callback(
#     dash.dependencies.Output('hover-data', 'children'),
#     [dash.dependencies.Input('indicator-graphic', 'hoverData')])
# def display_hover_data(hoverData):
#     return json.dumps(hoverData, indent=2)

@app.callback(
    dash.dependencies.Output('indicator-graphic', 'figure'),
    [dash.dependencies.Input('xaxis-column', 'value'),
     dash.dependencies.Input('yaxis-column', 'value'),
     dash.dependencies.Input('year--slider', 'value')])
def update_graph(xcolumn, ycolumn,
                 year_value):
    dff = df[df['Round'] == year_value]
    return scatter_plot(dff, xcolumn, ycolumn)

@app.callback(
    dash.dependencies.Output('chem_score', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value')])
def return_molecule_name(hoverData, year_value):
    dff = df[df['Round'] == year_value]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                curve_number = firstPoint['curveNumber']
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['category'] == category]
                score = dfff['Tanimoto Similarity Score'].iloc[point_number]
                tag = "Similarity score:  {:.2}".format(score)
                return tag

@app.callback(
    dash.dependencies.Output('chem_relative', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value')])
def return_molecule_name(hoverData, year_value):
    dff = df[df['Round'] == year_value]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                curve_number = firstPoint['curveNumber']
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['category'] == category]
                relative = dfff['Molecular Relative'].iloc[point_number]
                tag = "Molecular relative:  {}".format(relative)
                return tag

@app.callback(
    dash.dependencies.Output('chem_atoms', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value')])
def return_molecule_name(hoverData, year_value):
    dff = df[df['Round'] == year_value]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                curve_number = firstPoint['curveNumber']
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['category'] == category]
                atoms = dfff['Cation Heavy Atoms'].iloc[point_number]
                tag = "Heavy atoms:  {}".format(atoms)
                return tag
@app.callback(
    dash.dependencies.Output('chem_anion', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value')])
def return_molecule_name(hoverData, year_value):
    dff = df[df['Round'] == year_value]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                curve_number = firstPoint['curveNumber']
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['category'] == category]
                anion = dfff['Anion'].iloc[point_number]
                tag = "Anion:  {}".format(anion)
                return tag

@app.callback(
    dash.dependencies.Output('chem_smiles', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value')])
def return_molecule_name(hoverData, year_value):
    dff = df[df['Round'] == year_value]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                curve_number = firstPoint['curveNumber']
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['category'] == category]
                smiles = dfff['Salt Smiles'].iloc[point_number]
                tag = "SMILES:  {}".format(smiles.split(".")[0])
                return tag

@app.callback(
    dash.dependencies.Output('chem_img', 'src'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value')])
def display_image(hoverData, year_value):
    dff = df[df['Round'] == year_value]
    global encoded_image
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                curve_number = firstPoint['curveNumber']
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                #trace_name = app.layout['indicator-graphic'].figure['data']\
                #    [curve_number]['name']
                dfff = dff[dff['category'] == category]
                image_ID = \
                dfff['Salt Smiles'].iloc[
                    point_number].split(".")[0]
                image_filename= asset_directory + image_ID + '.png'
                encoded_image = base64.b64encode(
                    open(image_filename, 'rb').read())
                decoded_image = 'data:image/png;base64,{}'.format(encoded_image.decode())
                return decoded_image

external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "//fonts.googleapis.com/css?family=Raleway:400,300,600",
    "//fonts.googleapis.com/css?family=Dosis:Medium",
    "https://cdn.rawgit.com/plotly/dash-app-stylesheets/0e463810ed36927caf20372b6411690692f94819/dash-drug-discovery-demo-stylesheet.css"]

for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == '__main__':
    app.run_server()