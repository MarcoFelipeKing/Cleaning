#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 18:22:52 2020

@author: marcofking
"""


import pyabc as pyabc
from pyabc import (ABCSMC,
                   RV, Distribution,
                   MedianEpsilon,
                   LocalTransition)
from pyabc.visualization import plot_kde_2d, plot_data_callback
import matplotlib.pyplot as plt
import os
import tempfile
import numpy as np
#import scipy as sp
from scipy.integrate import odeint
import math
db_path = ("sqlite:///" +
           os.path.join(tempfile.gettempdir(), "test3.db"))


initial_contamination=59
measurement_data = np.array([19,5,5,2,9])
s=np.array([26,2.3,4.67,4.33,4.27])
precision=5000


measurement_times = np.array([1,2,4,8,24])#np.arange(len(measurement_data))  



def Distance(simulation, data):
    return np.absolute((data["Contamination"] - simulation["Contamination"])/data["sd"]).sum()

# def Distance(x,y,s):

#     # computes the Euclidean distance between two lists of the same length

#     if len(x) == len(y):

#         return math.sqrt(sum([(((x[i]-y[i])/s[i])**2) for i in range(len(x))]))

#     else:

#         return 100.00


def ode_model(contamination,t,r,C,d,g):
    Contamination = contamination;
    return(Contamination*r*(1-Contamination/C)-d*math.exp(-g*t)*Contamination)

# No Noise
def deterministic_run_NONOISE(parameters):#precision,initial_contamination,r,C,d,g):
    precision=5000
    tmax = 24
    time_space = np.linspace(0,tmax,precision+1)
    sim=odeint(ode_model,initial_contamination,time_space,args=(parameters["r"],parameters["C"],parameters["d"],parameters["g"]))
    #num_at_0=sim[int(precision*0.1/50.0)]
    num_at_1=sim[int(precision*1/tmax)]
    num_at_2=sim[int(precision*2/tmax)]
    num_at_4=sim[int(precision*4/tmax)]
    num_at_8=sim[int(precision*8/tmax)]
    num_at_24=sim[int(precision*24/tmax)]
    return{"Contamination":[num_at_1,num_at_2,num_at_4,num_at_8,num_at_24]}

def deterministic_run(parameters):#precision,initial_contamination,r,C,d,g):
    precision=5000
    tmax = 24
    time_space = np.linspace(0,tmax,precision+1)
    sim=odeint(ode_model,initial_contamination,time_space,args=(parameters["r"],parameters["C"],parameters["d"],parameters["g"]))
    #num_at_0=sim[int(precision*0.1/50.0)]
    num_at_1=sim[int(precision*1/tmax)]
    num_at_2=sim[int(precision*2/tmax)]
    num_at_4=sim[int(precision*4/tmax)]
    num_at_8=sim[int(precision*8/tmax)]
    num_at_24=sim[int(precision*24/tmax)]
    return{"Contamination":[num_at_1,num_at_2,num_at_4,num_at_8,num_at_24]+ sigma*np.random.randn(5)}

# def f(y, t0, theta1, theta2):
#     x1, x2 = y
#     dx1 = - theta1 * x1 + theta2 * x2
#     dx2 =   theta1 * x1 - theta2 * x2
#     return dx1, dx2
    
# def model(pars):
#     sol = sp.integrate.odeint(
#              f, init, measurement_times,
#              args=(pars["theta1"],pars["theta2"]))
#     return {"X_2": sol[:,1]}

# true_trajectory = model({"theta1": theta1_true,
#                          "theta2": theta2_true})["X_2"]

# plt.plot(true_trajectory, color="C0", label='Simulation')
# plt.scatter(measurement_times, measurement_data,
#             color="C1", label='Data')
# plt.xlabel('Time $t$')
# plt.ylabel('Measurement $Y$')
# plt.title('Conversion reaction: True parameters fit')
# plt.legend()
# plt.show()

# def distance(simulation, data):
#     return np.absolute(data["X_2"] - simulation["X_2"]).sum()

parameter_prior = Distribution(r=RV("uniform", 0.1, 4.0),
                               C=RV("uniform", 6.0, 10.0),
                               d=RV("uniform", 0.01, 4.0),
                               g=RV("uniform", 0.01, 4.0))

parameter_prior.get_parameter_names()



#Noisey model
# sigma=0.02
# acceptor = pyabc.StochasticAcceptor()
# kernel = pyabc.IndependentNormalKernel(var=sigma**2)
# eps = pyabc.Temperature()

# abc = pyabc.ABCSMC(deterministic_run, parameter_prior, kernel, eps=eps, acceptor=acceptor,population_size=100)
# abc.new(db_path,{"Contamination": measurement_data}) #This distance model assumes the name of the predicited and confirmed are the same
# history_acceptor = abc.run(max_nr_populations=10,minimum_epsilon=10)


#No Noise
abc = ABCSMC(models=deterministic_run_NONOISE,
              parameter_priors=parameter_prior,
              distance_function=Distance,
              population_size=50,
              transitions=LocalTransition(k_fraction=.5),
              eps=MedianEpsilon(500, median_multiplier=0.7))

abc.new(db_path, {"Contamination": measurement_data,"sd":s})
history = abc.run(minimum_epsilon=12, max_nr_populations=10)


from pyabc.visualization import plot_kde_matrix

df, w = history.get_distribution(m=0)
plot_kde_matrix(df, w);

# pyabc.visualization.plot_sample_numbers(history_acceptor)

# from pyabc.visualization import plot_kde_matrix

# df, w = history_acceptor.get_distribution(m=0)
# #plot_kde_matrix(df, w);
# df.hist(color='k', alpha=0.5, bins=25)

df.describe()

# See convergence of the distance function
pyabc.visualization.plot_distance_convergence(history)

# See convergence of the kernel width
fig, ax = plt.subplots()
for t in range(history.max_t + 1):
    df, w = history.get_distribution(m=0, t=t)
    pyabc.visualization.plot_kde_1d(
        df,
        w,
        xmin=0,
        xmax=1,
        x="lam",
        xname=r"$\lambda$",
        ax=ax,
        label=f"PDF t={t}",
    )
#Add a vertical line at lambda=0.49 +- 0.32 to 0.72 
#This is the experimental value predicted in King et al. 2020 with Kalanne
ax.axvline(0.49, color="k", linestyle="dashed")
ax.legend();

#Plot some trajectories

import pandas as pd
import operator

initial_contamination=59
measurement_data = np.array([59,19,5,5,2,9])
s=np.array([30,26,2.3,4.67,4.33,4.27])
precision=5000


measurement_times = np.array([0,1,2,4,8,24])

P=deterministic_run({"r":df["r"].mean(),"C":df["C"].mean(),"d":df["d"].mean(),"g":df["g"].mean()})
Pmin=deterministic_run({'r':df["r"].min(),"C":df["C"].min(),"d":df["d"].min(),"g":df["g"].min()})
Pmax=deterministic_run({'r':df["r"].max(),"C":df["C"].max(),"d":df["d"].max(),"g":df["g"].max()})

# create a vector of values between 0 and 5
x = np.linspace(0, 24, 6)

#Define new sd just for plotting to avoid SD value at 0


#Plot errobars of experimental data
plt.errorbar(measurement_times,measurement_data,yerr=s,fmt='o', color='black',label='Experimental data')

#Plot the model prediction
plt.plot(x,P["Contamination"],label="Model prediction",color='blue')

#Plot confidence intervals around the model prediction
plt.fill_between(x,Pmin["Contamination"],Pmax["Contamination"],alpha=0.2,color='blue')



#plt.fill_between(x, np.array(map(operator.sub, P["Contamination"], Pmin["Contamination"])), np.array(map(operator.add, P["Contamination"], Pmax["Contamination"])), color='b', alpha=.1)
plt.xlim(-1,6)
plt.ylabel("$\mu g$ recovered from finger \n after n contacts")
plt.xlabel("Number of contacts")
plt.legend(loc="upper left")


#save the plot
#plt.savefig("../Images/abc_prediction.png", dpi=600)

plt.show()


# fig = plt.figure(figsize=(10,8))
# for t in range(history_acceptor.max_t+1):
#     ax = fig.add_subplot(3, np.ceil(history_acceptor.max_t / 3), t+1)

#     ax = plot_kde_2d(
#         *history_acceptor.get_distribution(m=0, t=t), "r", "C","d","g",
#         xmin=0, xmax=15, numx=200, ymin=0, ymax=15, numy=200, ax=ax)
#     #ax.scatter([theta1_true], [theta2_true], color="C1",
#     #            label='$\Theta$ true = {:.3f}, {:.3f}'.format(
#     #                theta1_true, theta2_true))
#     ax.set_title("Posterior t={}".format(t))

#     ax.legend()
# fig.tight_layout()



# _, ax = plt.subplots()
# for t in range(history_acceptor.max_t + 1):
#     pyabc.visualization.plot_kde_1d_highlevel(
#         history_acceptor, x="Contamination", t=t,
#         refval=measurement_data, refval_color='grey',
#         xmin=0, xmax=15, ax=ax, numx=50, label=f"Iteration {t}")
# ax.legend()
# plt.show()


# # History of acceptances
# pyabc.visualization.plot_sample_numbers(history_acceptor, labels="noisy")
# plt.show()

# ####

# _, ax = plt.subplots()

# def plot_data(sum_stat, weight, ax, **kwargs):
#     """Plot a single trajectory"""
#     ax.plot(measurement_times, sum_stat['Contamination'], color='grey', alpha=0.1)
    
# def plot_mean(sum_stats, weights, ax, **kwargs):
#     """Plot mean over all samples"""
#     weights = np.array(weights)
#     weights /= weights.sum()
#     data = np.array([sum_stat['Contamination'] for sum_stat in sum_stats])
#     mean = (data * weights.reshape((-1, 1))).sum(axis=0)
#     ax.plot(measurement_times, mean, color='C2', label='Sample mean')
    
# ax = plot_data_callback(h, plot_data, plot_mean, ax=ax)

# #plt.plot(true_trajectory, color="C0", label='Simulation')
# plt.errorbar([0,1,2,4,8,24],np.append(initial_contamination,measurement_data), yerr=[92,26,2.3,4.67,4.33,4.2], fmt='x',color="Teal")
# plt.xlabel('Time $t$')
# plt.ylabel('CFU ')
# plt.title('Time(h) after cleaning')
# plt.legend()
# plt.show()