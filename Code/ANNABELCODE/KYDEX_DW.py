#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Built-in libraries
import os
import sys
from io import StringIO

# Third-party libraries
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from tqdm import tqdm
import warnings
import seaborn as sns
import matplotlib.pyplot as plt

from contextlib import redirect_stderr


# Constants and Configuration
EPSILON = 250
NUM_SAMPLES = 1000
DATA_PATH = '../../Annabel_data/FINAL_cleaning_data.csv'

def load_and_process_data():
    """Load the data, filter it, and group by time."""
    df = pd.read_csv(DATA_PATH)
    df = df[(df['Surface'] == "KYDEX") & (df['Cleaning'] == "DW") & (df['Time'] > -.5)]
    grouped_data = df.groupby('Time').agg({'Count': ['mean', 'std']}).reset_index()
    grouped_data.columns = ['Time', 'mean', 'std']
    return grouped_data

def logistic_decay_model(params, times, grouped_data):
    """... [Your existing docstring or a new one describing the function] ..."""

    np.seterr(invalid='ignore')
    warnings.filterwarnings('ignore')
    #surface_key = surface.replace(' ', '_')  # Replace space with underscore
    r = params['r']#params[surface_key + '_r_' + cleaning]
    K = params['K']#params[surface_key + '_K_' + cleaning]
    m = params['m']#params[surface_key + '_m_' + cleaning]
    g = params['g']#params[surface_key + '_g_' + cleaning]
    C0 = params['C0']#params[surface_key + '_C0' ]  # Use the initial condition from the parameters
    times = np.array([0,1,2,4])
    def dCdt(C, t):
        cleaning_effect = m * np.exp(-g * t) * C if t >= 0 else 0
        return r * C * (1 - C / K) - cleaning_effect

    # Suppress the specified warnings
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    warnings.filterwarnings('ignore', message='t + h = t')
    warnings.filterwarnings('ignore', message='lsoda--  warning..internal t')
    C = odeint(dCdt, C0, times)
    return C.flatten()

def sample_params():
    """Sample parameters for the model."""
    r = 0 + 25 * np.random.rand(1)
    K = 0 + 3E4 * np.random.rand(1)
    m = 0 + 20 * np.random.rand(1)
    g = 0 + 25 * np.random.rand(1)
    C0 = 500 + 2.5E3 * np.random.randn(1)
    return r, K, m, g, C0


def main():
    grouped_data = load_and_process_data()
    
    accepted_params = []
    distances = []  # List to store the distances
    
    with tqdm(total=NUM_SAMPLES, desc="Sampling", ncols=100) as pbar:
        # Suppress stderr warnings
        with open(os.devnull, 'w') as fnull:
            with redirect_stderr(fnull):
                while len(accepted_params) < NUM_SAMPLES:
                    r, K, m, g, C0 = sample_params()
                    params = {'r': r, 'K': K, 'm': m, 'g': g, 'C0': C0}
                    sim = logistic_decay_model(params, grouped_data['Time'].values, grouped_data)
                    dist = np.sum(abs(sim - grouped_data['mean']))
                    
                    if dist < EPSILON:
                        accepted_params.append(params)
                        distances.append(dist)
                        pbar.update(1)

    # Convert the accepted parameters and distances to a DataFrame and save to CSV
    df_accepted = pd.DataFrame(accepted_params)
    df_accepted['distance'] = distances
    df_accepted.to_csv('accepted_parameters.csv', index=False)

if __name__ == "__main__":
    main()

def plot_and_save_pairplot():
    # Load the accepted parameters from CSV
    df_accepted = pd.read_csv('accepted_parameters.csv')
    
    # Drop the distance column for the pairplot
    df_plot = df_accepted.drop(columns=['distance'])
    
    # Explicitly specify the columns in pairplot
    columns_to_plot = ['r', 'K', 'm', 'g', 'C0']
    sns.pairplot(df_plot[columns_to_plot])
    plt.savefig('pairplot.png')

# Call the function after main
plot_and_save_pairplot()

def plot_best_prediction():
    # Load the accepted parameters from CSV
    df_accepted = pd.read_csv('accepted_parameters.csv')
    
    # Reset the index and convert list-type values to scalar values
    df_accepted.reset_index(inplace=True)
    for col in df_accepted.columns:
        df_accepted[col] = df_accepted[col].apply(lambda x: eval(x)[0] if isinstance(x, str) and '[' in x else x)
    
    # Sort by distance and get the best parameters
    df_accepted.sort_values(by='distance', inplace=True)
    best_params = df_accepted.iloc[0].drop('distance').to_dict()
    
    # Use the best parameters to generate predictions
    times = np.array([0,1,2,4])
    predictions = logistic_decay_model(best_params, times, grouped_data)
    
    # Plotting the results
    plt.figure(figsize=(10,6))
    
    # Plot experimental data with error bars
    plt.errorbar(grouped_data['Time'].values, grouped_data['mean'].values, yerr=grouped_data['std'].values, 
                 fmt='o', color='black', label='Experimental data')
    
    # Plot predictions
    plt.plot(times, predictions, '-r', label='Best prediction')
    
    plt.xlabel('Time [h]]')
    plt.ylabel('CFUs')
    plt.title('Best Model Prediction vs Experimental Data')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.savefig('best_prediction.png')
    plt.show()

# Load and process your data
grouped_data = load_and_process_data()

# Call the plotting function
plot_best_prediction()