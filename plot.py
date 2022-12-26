import numpy as np
import plotly.express as px
import pandas as pd
import sys

def read_file(file_name):
    f = open(file_name)
    res_x = []
    res_y = []
    res_z = []
    value = []
    size = []
    for line_number, line in enumerate(f):
        words = line.strip().split()
    
        res_x.append(float(words[0]))
        res_y.append(float(words[1]))
        res_z.append(float(words[2]))
        value.append(float(words[3]))
        size.append(0.01)
    data = {"x":res_x, "y":res_y, "z":res_z, "value":value, "size":size}
    return data

def print_image(file_name):
    data = read_file(file_name)
    df = pd.DataFrame(data)

    fig = px.scatter_3d(df, x='x', y='y', z='z',
                    color='value', size="size", opacity=0.8)
    fig.write_image(file_name+".pdf", width=2000, height=2000, scale=3)

if __name__ == "__main__":
    print_image("an_sol_0_0_0")
    print_image("calc_sol_0_0_0")
    print_image("calc_sol_err_0_0_0")
