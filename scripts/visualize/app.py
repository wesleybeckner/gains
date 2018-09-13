import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import numpy as np
import os
import base64

app = dash.Dash()
asset_directory = os.getcwd() + '/assets/'
df = pd.read_csv("assets/vizapp.csv")

# create controls

dataset_options = [{'label': 'first series', 'value': 1},
                   {'label': 'second series', 'value': 2},
                   {'label': 'third series', 'value': 3}]
starting_series = 1
starting_round = 1

dff = df.loc[df['Round'] == starting_round]
cols = []
for item in df.columns:
    if "density" in item:
        cols.append(item)
    elif "cpt" in item:
        cols.append(item)
    elif "Score" in item:
        cols.append(item)
available_indicators = cols
xcolumn = available_indicators[0]
ycolumn = available_indicators[1]

pareto_xs = [730.0,
             1635.0,
             310.0,
             208.0,
             502.0,
             1635.0,
             502.0,
             298.0,
             232.0,
             208.0,
             232.0,
             298.0,
             522.0,
             730.0,
             522.0,
             310.0]

pareto_ys = [1551.0,
             887.0,
             1470.0,
             1388.0,
             1014.0,
             887.0,
             1014.0,
             1047.0,
             1153.0,
             1388.0,
             1153.0,
             1047.0,
             1536.0,
             1551.0,
             1536.0,
             1470.0]


def scatter_plot(
        pareto_xs,
        pareto_ys,
        dff,
        xcolumn,
        ycolumn,
        size=dff["Cation Heavy Atoms"],
        plot_type='Scatter'):
    trace0 = [dict(
        x=dff[dff['Category'] == i][xcolumn],
        y=dff[dff['Category'] == i][ycolumn],
        xlabel=xcolumn,
        ylabel=ycolumn,
        customdata=dff[dff['Category'] == i]['Category'],
        smiles=dff['Salt Smiles'],
        score=dff[dff['Category'] == i]['Tanimoto Similarity Score'],
        atoms=dff['Cation Heavy Atoms'],
        relative=dff['Molecular Relative'],
        anion=dff['Anion'],
        category=dff['Category'],
        name=i,
        mode='markers',
        marker=dict(
            sizeref=45,
            sizemode='diameter',
            opacity=0.5,
            size=size * 60,
        ),
        type=plot_type,
    ) for i in dff.Category.unique()
    ]
    trace_pareto = [dict(
        x=pareto_xs,
        y=pareto_ys,
        mode='Scatter',
        dash='dash',
        connectgaps=True,
        fillcolor='#FAEBFC',
        line=dict(
            color='rgba(32, 32, 32, .6)',
            width=1
        ),
        name='pareto front',
        hoverinfo='skip',
        type='Scatter',
        visible='legendonly'
    )
    ]
    layout = dict(
        margin=dict(r=15, t=15),  # , l=40, b=40),
        xaxis=dict(
            title=xcolumn,
        ),
        yaxis=dict(
            title=ycolumn,
        ),
        showlegend=True,
    )

    return dict(data=trace0 + trace_pareto, layout=layout)


# for side panel

smiles = df.loc[df['Round'] == starting_round]['Salt Smiles'].iloc[1]
score = df.loc[df['Round'] == starting_round]['Tanimoto Similarity Score'].iloc[1]
atoms = df.loc[df['Round'] == starting_round]['Cation Heavy Atoms'].iloc[1]
relative = df.loc[df['Round'] == starting_round]['Molecular Relative'].iloc[1]
anion = df.loc[df['Round'] == starting_round]['Anion'].iloc[1]
category = df.loc[df['Round'] == starting_round]['Category'].iloc[1]
image_ID = df.loc[df['Round'] == starting_round]['Salt Smiles'].iloc[1].split(".")[0]

logo_filename = asset_directory + 'gains.png'
image_filename = asset_directory + image_ID + '.png'
encoded_logo = base64.b64encode(open(logo_filename, 'rb').read())
encoded_image = base64.b64encode(open(image_filename, 'rb').read())
FIGURE = scatter_plot(pareto_xs, pareto_ys, dff, xcolumn, ycolumn)

app.layout = html.Div([
    # Row 1: Header and Intro text
    html.Div([

        html.Img(src='data:image/png;base64,{}'.format(encoded_logo.decode()),
                 style={
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
                    'HOVER over a salt in the graph to the right to see its '
                    'structure and metadata to the left'),
                html.P(
                    'SELECT parameters in the dropdown '
                    'menus to set axes and populate the graph')
#                dcc.Checklist(
#                    id='lock_selector',
#                    options=[
#                        {'label': 'Lock camera', 'value': 'locked'}
#                    ],
#                    values=[],
#                )
            ], style={'margin-left': '10px'}),
        ], className='twelve columns')

    ], className='row'),

    # Row 2: Hover Panel and Graph
    html.Div([
        html.Div([
           # dcc.RadioItems(
           #     id='series_selector',
           #     options=[
           #         {'label': 'All ', 'value': 'all'},
           #         {'label': 'None ', 'value': 'none'}
           #     ],
           #     value='all',
           #     labelStyle={'display:': 'inline-block'}
           # ),
            dcc.Dropdown(
                id='search_series',
                options=dataset_options,
                multi=True,
                value=['first series']
            ),
            html.Img(id='chem_img', src='data:image/png;base64,{}'.format(encoded_image.decode()),
                     style={
                         'height': 'auto',
                         'width': '100%',
                         'float': 'top',
                         'position': 'relative',
                         'bottom': '0px',
                     }),
            html.Br(),

            html.P(
                'Description',
                style={
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

        ], className='three columns',
            style={'margin-left': '10px'}),
        # ], className='three columns', style=dict(height='300px', left='10px')),
        html.Div([
            dcc.Dropdown(
                id='xaxis-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='predicted cpt'
            ),

        ], className='four columns'),

        html.Div([
            dcc.Dropdown(
                id='yaxis-column',
                options=[{'label': i, 'value': i} for i in available_indicators],
                value='predicted density'
            ),

        ], className='four columns'),
    ]),
    html.Div([
        dcc.Graph(id='indicator-graphic',
                  figure=FIGURE),

        dcc.RangeSlider(
            id='year--slider',
            min=df['Round'].min(),
            max=df['Round'].max(),
            value=[df['Round'].min(), df['Round'].max()],
            step=None,
            marks={str(year): str(year) for year in df['Round'].unique()}
        ),
        html.P(' '),
        html.P('Selection mean error: {:.3}'.format(score),
               id='selection-summary'),
    ], className='eight columns'
    ),
])

@app.callback(
    dash.dependencies.Output('selection-summary', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'selectedData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def return_molecule_name(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    scores = []
    if hoverData is not None:
        if 'points' in hoverData:
            for datapoint in hoverData['points']:
                point_number = datapoint['pointNumber']
                category = datapoint['customdata']
                dfff = dff[dff['Category'] == category]
                score = dfff['Tanimoto Similarity Score'].iloc[point_number]
                scores.append(score)
            tag = "Average similarity score:  {:.2}".format(np.average(scores))
            return tag


@app.callback(
    dash.dependencies.Output('indicator-graphic', 'figure'),
    [dash.dependencies.Input('xaxis-column', 'value'),
     dash.dependencies.Input('yaxis-column', 'value'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def update_graph(xcolumn, ycolumn,
                 year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    return scatter_plot(pareto_xs, pareto_ys, dff, xcolumn, ycolumn)


@app.callback(
    dash.dependencies.Output('chem_score', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def return_molecule_name(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['Category'] == category]
                score = dfff['Tanimoto Similarity Score'].iloc[point_number]
                tag = "Similarity score:  {:.2}".format(score)
                return tag


@app.callback(
    dash.dependencies.Output('chem_relative', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def return_molecule_name(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['Category'] == category]
                relative = dfff['Molecular Relative'].iloc[point_number]
                tag = "Molecular relative:  {}".format(relative)
                return tag


@app.callback(
    dash.dependencies.Output('chem_atoms', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def return_molecule_name(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['Category'] == category]
                atoms = dfff['Cation Heavy Atoms'].iloc[point_number]
                tag = "Heavy atoms:  {}".format(atoms)
                return tag


@app.callback(
    dash.dependencies.Output('chem_anion', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def return_molecule_name(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['Category'] == category]
                anion = dfff['Anion'].iloc[point_number]
                tag = "Anion:  {}".format(anion)
                return tag


@app.callback(
    dash.dependencies.Output('chem_smiles', 'children'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def return_molecule_name(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['Category'] == category]
                smiles = dfff['Salt Smiles'].iloc[point_number]
                tag = "SMILES:  {}".format(smiles.split(".")[0])
                return tag


@app.callback(
    dash.dependencies.Output('chem_img', 'src'),
    [dash.dependencies.Input('indicator-graphic', 'hoverData'),
     dash.dependencies.Input('year--slider', 'value'),
     dash.dependencies.Input('search_series', 'value')])
def display_image(hoverData, year_value, series_value):
    dff = df[df['Series'].isin(series_value)
             & df['Round'].between(year_value[0], year_value[1])]
    global encoded_image
    if hoverData is not None:
        if 'points' in hoverData:
            firstPoint = hoverData['points'][0]
            if 'pointNumber' in firstPoint:
                point_number = firstPoint['pointNumber']
                category = firstPoint['customdata']
                dfff = dff[dff['Category'] == category]
                image_ID = \
                    dfff['Salt Smiles'].iloc[
                        point_number].split(".")[0]
                image_filename = asset_directory + image_ID + '.png'
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
