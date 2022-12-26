import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import sys

val_mpi_1 = {"128":[1, 3.73, 5.15, 5.85, 12.10], "256":[1, 3.94, 6.53, 11.81, 20.85], 
             "512":[1, 3.96, 7.39, 13.93, 27.04], "MPI cpus num":[1, 4, 8, 16, 32]}

val_mpi_2 = {"128":[1,3.96,6.57,6.44,8.07], "256":[1, 3.98, 6.76, 13.00, 22.54], 
             "512":[1, 3.96, 7.77, 14.59, 28.37], "MPI cpus num":[1, 4, 8, 16, 32]}

val_omp_1 = {"128":[1, 6.20, 10.72, 17.67], "256":[1, 5.30, 10.36, 19.57], 
             "512":[1, 2.27, 4.42, 8.55], "MPI cpus num":[1, 2, 4, 8]}

val_omp_2 = {"128":[1, 6.61, 11.61, 19.37], "256":[1, 5.03, 9.89, 18.67], 
             "512":[1, 1.94, 3.87, 7.64], "MPI cpus num":[1, 2, 4, 8]}

def print_image(data, file_name):
    df = pd.DataFrame(data)

    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=df["MPI cpus num"], y=df["128"],
                    mode='lines+markers',
                    name="сетка 128*128*128"))
    fig.add_trace(go.Scatter(x=df["MPI cpus num"], y=df["256"],
                    mode='lines+markers',
                    name="сетка 256*256*256"))
    fig.add_trace(go.Scatter(x=df["MPI cpus num"], y=df["512"],
                    mode='lines+markers',
                    name="сетка 512*512*512"))
    fig.write_image(file_name+".pdf", width=500, height=500, scale=5)

if __name__ == "__main__":
    print_image(val_mpi_1,"mpi_l_0")
    print_image(val_mpi_2,"mpi_l_1")
    print_image(val_omp_1,"omp_l_0")
    print_image(val_omp_2,"omp_l_1")
