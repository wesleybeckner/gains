import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import re
from math import log
import numpy as np

app = dash.Dash()
md_plot = False
if md_plot:
    log = pd.read_csv("../../salty/scripts/salt_log.csv")
    predictions = []
    calculations = []
    errors = []
    for j in range(log.shape[0]):
        prediction = log["Model Prediction"][j]
        calculation = log["MD Calculation"][j]
        error = log["Error"][j]
        predictions.append(re.findall("\d+\.\d+", prediction))
        calculations.append(re.findall("\d+\.\d+", calculation))
        errors.append(re.findall("\d+\.\d+", error))
    predictions = pd.DataFrame(predictions)
    predictions.rename(columns=lambda x: "prediction {}".format(x), inplace=True)
    calculations = pd.DataFrame(calculations)
    calculations.rename(columns=lambda x: "calculation {}".format(x), inplace=True)
    errors = pd.DataFrame(errors)
    errors.rename(columns=lambda x: "error {}".format(x), inplace=True)
    new = pd.DataFrame(log[log.columns[:-3]])
    df = pd.concat([new,predictions,calculations,errors], axis=1)

df = pd.read_csv("vizapp.csv")


app.layout = html.Div([
    dcc.Graph(
        id='cpt-vs-density',
        figure={
            'data': [
                go.Scatter(
                    x=df[df['category'] == i]['Heat capacity at constant pressure, J/K/mol_mean'],
                    y=df[df['category'] == i]['Specific density, kg/m<SUP>3</SUP>_mean'],
                    text=df[df['category'] == i]['smiles-anion'],
                    mode='markers',
                    opacity=0.7,
                    marker={
                        'size': 15+3*np.log(df[df['category'] == i]['count']),
                        'line': {'width': 0.5, 'color': 'white'}
                    },
                    name=i
                ) for i in df.category.unique()
            ],
            'layout': go.Layout(
                xaxis={'title': 'Heat Capacity'},
                yaxis={'title': 'Density'},
                margin={'l': 40, 'b': 40, 't': 10, 'r': 10},

                #legend={'x': 0, 'y': 1},
                hovermode='closest'
            )
        }
    )
])

if __name__ == '__main__':
    app.run_server()
