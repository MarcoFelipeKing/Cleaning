import scipy.optimize as opt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

# Import data# import data using pandas


data = pd.read_csv('Annabel_data/FINAL_cleaning_data_NO_NIL.csv')

# 1. Impute missing values in 'Count' column with 10,000
data['Count'].fillna(10000, inplace=True)

# filter data t>=0
#data = data[data['Time'] >= 0]
def subtract_control(data, control):
    """
    Subtract control data from the given dataset.
    
    Parameters:
    - data: DataFrame containing the data to be adjusted.
    - control: DataFrame containing the control data.
    
    Returns:
    - Adjusted DataFrame.
    """
    data_copy = data.copy()  # Creating an explicit copy to avoid SettingWithCopyWarning
    
    # Find common times between data and control
    common_times = set(data['Time']).intersection(set(control['Time']))
    
    # Subtracting control mean and adjusting standard deviation for common times
    for time in common_times:
        data_index = data_copy['Time'] == time
        control_values = control[control['Time'] == time]
        
        data_copy.loc[data_index, 'Count mean'] -= control_values['Count mean'].values
        data_copy.loc[data_index, 'Count std'] = np.sqrt(
            data_copy.loc[data_index, 'Count std'].values**2 + control_values['Count std'].values**2
        )
    
    return data_copy


# 1. Group the data to calculate mean and standard deviation
grouped = data.groupby(['Time', 'Surface', 'Cleaning']).agg(
    {'Count':['mean', 'std']}).reset_index()

# Flatten the multi-level column index
grouped.columns = [' '.join(col).strip() for col in grouped.columns.values]

# 2. Rename the columns for easier access
grouped.rename(columns={'Count mean': 'Count mean', 'Count std': 'Count std'}, inplace=True)


# 3. Split the data based on surface and cleaning method
grouped_data = {}
for surface in ['KYDEX', 'Stainless Steel']:
    for cleaning in ['CON', 'ALC', 'DW']:
        key = f"{surface}_{cleaning}"
        grouped_data[key] = grouped[(grouped['Surface'] == surface) & (grouped['Cleaning'] == cleaning)]

# Subtract the control data from each cleaning method for each material
grouped_data_adjusted = {}

for key in grouped_data:
    if 'CON' not in key:  # We don't adjust the control data itself
        material, cleaning_method = key.split("_")
        control_key = f"{material}_CON"
        grouped_data_adjusted[key] = subtract_control(grouped_data[key], grouped_data[control_key])

# Adding back the control data to the dictionary
for key in grouped_data:
    if 'CON' in key:
        grouped_data_adjusted[key] = grouped_data[key]

# Display the processed data

# Define the logistic growth decay function
def logistic_growth_decay(C, t, r, K, m, g):
    cleaning_effect = m * np.exp(-g * t) if t >= 0 else 0
    dCdt = r * C * (1 - C / K) - cleaning_effect * C
    return dCdt

# Function to fit the logistic growth decay model to data
def fit_logistic_growth_decay(data):
    # Objective function to minimize
    def residuals(params):
        r, K, m, g = params
        simulated = odeint(logistic_growth_decay, data['Count mean'].iloc[0], data['Time'], args=(r, K, m, g)).flatten()
        return simulated - data['Count mean']
    
    # Initial guess for parameters
    initial_guess = [0.1, 20000, 0.5, 0.1]
    result = least_squares(residuals, initial_guess)
    return result.x

# Function to simulate and plot results
def simulate_and_plot(data, params, title, ax):
    times = np.linspace(data['Time'].min(), data['Time'].max(), 100)
    result = odeint(logistic_growth_decay, data['Count mean'].iloc[0], times, args=tuple(params))
    
    # Plotting
    ax.errorbar(data['Time'].values, data['Count mean'].values, yerr=data['Count std'].values, fmt='o', label='Experimental Data')
    ax.plot(times, result, label='Fitted Logistic Growth Decay')
    
    # Adding parameter annotations
    ax.text(0.1, 0.95, f"r = {params[0]:.4f}", transform=ax.transAxes)
    ax.text(0.1, 0.90, f"K = {params[1]:.4f}", transform=ax.transAxes)
    ax.text(0.1, 0.85, f"m = {params[2]:.4f}", transform=ax.transAxes)
    ax.text(0.1, 0.80, f"g = {params[3]:.4f}", transform=ax.transAxes)

    ax.set_title(title)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Bacterial Count')
    ax.legend()
    ax.grid(True)


def fit_and_plot_surface_cleaning_combinations(grouped_data, combinations, subtract_control_flag=True, time_range=None):
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    axes = axs.ravel()
    fitted_params = {}
    
    for i, (title, (surface, cleaning)) in enumerate(combinations.items()):
        key = f"{surface}_{cleaning}"
        
        # Filter data based on time range
        if time_range:
            start_time, end_time = time_range
            data_subset = grouped_data[key][(grouped_data[key]['Time'] >= start_time) & (grouped_data[key]['Time'] <= end_time)]
        else:
            data_subset = grouped_data[key]
        
        # Subtract control if the flag is set
        if subtract_control_flag and 'CON' not in key:
            control_key = f"{surface}_CON"
            if control_key in grouped_data and not data_subset['Time'].lt(0).any():  # Ensure control data exists and time is not < 0
                data_subset = subtract_control(data_subset, grouped_data[control_key])
        
        params = fit_logistic_growth_decay(data_subset)
        fitted_params[title] = params
        simulate_and_plot(data_subset, params, title, axes[i])

    plt.tight_layout()
    plt.show()

    return fitted_params





# Define combinations
combinations = {
    'KYDEX: ALC': ('KYDEX', 'ALC'),
    'KYDEX: DW': ('KYDEX', 'DW'),
    'Stainless Steel: ALC': ('Stainless Steel', 'ALC'),
    'Stainless Steel: DW': ('Stainless Steel', 'DW')
}

# Fit the model and plot results
fitted_parameters = fit_and_plot_surface_cleaning_combinations(grouped_data, combinations, subtract_control_flag=False, time_range=(0, 8))
#fit_and_plot_surface_cleaning_combinations(grouped_data_adjusted, combinations)

# Display the results
for combination, params in fitted_parameters.items():
    print(f"{combination} Parameters:")
    print(f"r = {params[0]:.4f}")
    print(f"K = {params[1]:.4f}")
    print(f"m = {params[2]:.4f}")
    print(f"g = {params[3]:.4f}")
    print("-" * 40)
