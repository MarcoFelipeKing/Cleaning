import numpy as np

from scipy.integrate import odeint

from scipy.stats import median_abs_deviation 

from tqdm import tqdm
import pandas as pd
import os
import warnings
import scipy.integrate

import os
import sys

import numpy as np
import scipy.integrate as integrate
from numpy import pi
import contextlib

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextlib.contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    """
    https://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

import contextlib




def logistic_decay_model(params, times, surface, cleaning):
    np.seterr(invalid='ignore')
    warnings.filterwarnings('ignore')
    surface_key = surface.replace(' ', '_')  # Replace space with underscore
    r = params[surface_key + '_r_' + cleaning]
    K = params[surface_key + '_K_' + cleaning]
    m = params[surface_key + '_m_' + cleaning]
    g = params[surface_key + '_g_' + cleaning]
    C0 = params[surface_key + '_C0' ]  # Use the initial condition from the parameters

    def dCdt(C, t):
        cleaning_effect = m * np.exp(-g * t) * C if t >= 0 else 0
        return r * C * (1 - C / K) - cleaning_effect

    # Suppress the specified warnings
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    warnings.filterwarnings('ignore', message='t + h = t')
    warnings.filterwarnings('ignore', message='lsoda--  warning..internal t')
    C = odeint(dCdt, C0, times)
    return C.flatten()


def richardson_model(params, times, surface, cleaning):
    np.seterr(invalid='ignore')
    surface_key = surface.replace(' ', '_')
    r = params[surface_key + '_r']
    alpha = params[surface_key + '_alpha_' + cleaning]
    C0 = params[surface_key + '_C0_' + cleaning]

    def dCdt(C, t):
        return r * C**alpha

    C = odeint(dCdt, C0, times) #+ np.random.normal(0, 1, len(times))
    return C.flatten()


def gompertz_model(params, times, surface, cleaning):
    np.seterr(invalid='ignore')
    surface_key = surface.replace(' ', '_')
    r = params[surface_key + '_r']
    K = params[surface_key + '_K']
    C0 = params[surface_key + '_C0_' + cleaning]

    def dCdt(C, t):
        return r * np.log(K / C) * C 

    C = odeint(dCdt, C0, times) #+ np.random.normal(0, 1, len(times))
    return C.flatten()



def quorum_sensing_model(params, times, surface, cleaning):
    # Suppressing the warning
    np.seterr(invalid='ignore')
    surface_key = surface.replace(' ', '_')  # Replace space with underscore
    r = params[surface_key + '_r']
    K = params[surface_key + '_K']
    m = params[surface_key + '_m_' + cleaning]
    g = params[surface_key + '_g_' + cleaning]
    C0 = params[surface_key + '_C0_' + cleaning]  # Use the initial condition from the parameters
    q = params[surface_key + '_q_' + cleaning]
    Cq = params[surface_key + '_Cq_' + cleaning]
    n = params[surface_key + '_n_' + cleaning]

    def dCdt(C, t):
        quorum_sensing_term = q * (C ** n) / ((C ** n) + (Cq ** n))
        return r * C * (1 - C / K) - m * np.exp(-g * t) * C + quorum_sensing_term 

    C = odeint(dCdt, C0, times)
    return C.flatten() + np.random.normal(0, 10, len(times))

def quorum_sensing_autoinducer_model(params, times, surface, cleaning):
    # Suppressing the warning
    np.seterr(invalid='ignore')
    surface_key = surface.replace(' ', '_')  # Replace space with underscore
    r = params[surface_key + '_r']
    K = params[surface_key + '_K']
    C0 = params[surface_key + '_C0_' + cleaning]

    if cleaning == 'CON':
        # For CON cleaning, only fit the logistic growth model
        def logistic_growth(C, t):
            return r * C * (1 - C / K)
        C = odeint(logistic_growth, C0, times)
        return C.flatten()

    # For other cleaning types, include quorum sensing and cleaning terms
    m = params[surface_key + '_m_' + cleaning]
    g = params[surface_key + '_g_' + cleaning]
    p = params[surface_key + '_p_' + cleaning]
    d = params[surface_key + '_d_' + cleaning]
    n = params[surface_key + '_n_' + cleaning]
    K_QS = params[surface_key + '_K_QS_' + cleaning]
    A0 = params[surface_key + '_A0_' + cleaning]  # Initial condition for autoinducer

    def system(y, t):
        C, A = y
        quorum_sensing_term = 0 if t < 0 else A**n / (A**n + K_QS**n)
        #print('quorum_sensing_term: ', quorum_sensing_term)
        cleaning_term = 0 if t < 0 else m * np.exp(-g * t) * C
        #print('cleaning_term: ', cleaning_term)
        logistic_model = r * C * (1 - C / K)
        #print('logistic_model: ', logistic_model)
        dCdt = quorum_sensing_term * logistic_model - cleaning_term
        dAdt = p * C - d * A
        return [dCdt, dAdt]

    result = odeint(system, [C0, A0], times)
    C = result[:, 0]
    # if any NANs then replace with 10^6
    #C = np.where(np.isnan(C), 10**6, C)
    return C.flatten() + np.random.normal(0, 5, len(times))


# PRIORS: 

def sample_from_prior():
    priors = {
        
        'KYDEX_q_ALC': np.random.uniform(0, 12),
        'KYDEX_q_CON': 0,#np.random.uniform(0, 200),
        'KYDEX_q_DW': np.random.uniform(0, 11),
        'KYDEX_n_ALC': np.random.uniform(0, 5),
        'KYDEX_n_CON': 0,#np.random.uniform(0, 200),
        'KYDEX_n_DW': np.random.uniform(0, 5),
        'KYDEX_Cq_ALC': np.random.uniform(0, 12000),
        'KYDEX_Cq_CON': 0,#np.random.uniform(0, 200),
        'KYDEX_Cq_DW': np.random.uniform(0, 12000),
        # {'KYDEX': {'C0': 0.545593923608381, 'r': -0.20455899876837547, 'K': 12723834423.153332}, 
        # 'Stainless Steel': {'C0': 0.08468304373538209, 'r': -0.23263366414248537, 'K': 10180567735.662352}}
        

        'KYDEX_r': np.random.uniform(-2, 100),
        'KYDEX_K': np.random.uniform(10, 3E4),

        'KYDEX_m_ALC': np.random.uniform(0, 100),
        'KYDEX_g_ALC': np.random.uniform(-20, 30),
        'KYDEX_m_CON': 0,#np.random.uniform(0, 0),
        'KYDEX_g_CON': 0,#np.random.uniform(0, 0),
        'KYDEX_m_DW': np.random.uniform(0, 100),
        'KYDEX_g_DW': np.random.uniform(-20, 30),

        
        'KYDEX_alpha_DW': np.random.uniform(0, 20),
        'KYDEX_alpha_ALC': np.random.uniform(0, 20),
        'KYDEX_alpha_CON': np.random.uniform(0, 20),

        'Stainless_Steel_alpha_CON': 0,
        'Stainless_Steel_q_ALC': np.random.uniform(0, 10),
        'Stainless_Steel_q_CON': 0,#np.random.uniform(0, 200),
        'Stainless_Steel_q_DW': np.random.uniform(0, 10),
        'Stainless_Steel_n_ALC': np.random.uniform(0, 5),
        'Stainless_Steel_n_CON': 0,#np.random.uniform(0, 200),
        'Stainless_Steel_n_DW': np.random.uniform(0, 5),
        'Stainless_Steel_Cq_ALC': np.random.uniform(0, 1000),
        'Stainless_Steel_Cq_CON': 0,#np.random.uniform(0, 200),
        'Stainless_Steel_Cq_DW': np.random.uniform(0, 1000),
        
        'Stainless_Steel_r': np.random.uniform(-2, 100),
        'Stainless_Steel_K': np.random.uniform(10, 3E4),
        'Stainless_Steel_m_ALC': np.random.uniform(10, 100),
        'Stainless_Steel_g_ALC': np.random.uniform(-20, 30),
        'Stainless_Steel_m_CON': 0,#np.random.uniform(0, 0),
        'Stainless_Steel_g_CON': 0,#np.random.uniform(0, 0),
        'Stainless_Steel_m_DW': np.random.uniform(0, 100),
        'Stainless_Steel_g_DW': np.random.uniform(-20, 30),

        
        'Stainless_Steel_alpha_DW': np.random.uniform(0, 200),
        'Stainless_Steel_alpha_ALC': np.random.uniform(0, 200),

                # New priors for autoinducer dynamics
        'KYDEX_p_ALC': np.random.uniform(0, 10),
        'KYDEX_p_CON': np.random.uniform(0, 10),
        'KYDEX_p_DW': np.random.uniform(0, 10),
        'KYDEX_d_ALC': np.random.uniform(0, 10),
        'KYDEX_d_CON': np.random.uniform(0, 10),
        'KYDEX_d_DW': np.random.uniform(0, 1),
        'KYDEX_K_QS_ALC': np.random.uniform(0, 1000),
        'KYDEX_K_QS_CON': np.random.uniform(0, 1000),
        'KYDEX_K_QS_DW': np.random.uniform(0, 1000),
        'KYDEX_A0_ALC': np.random.uniform(0, 100),
        'KYDEX_A0_CON': np.random.uniform(0, 100),
        'KYDEX_A0_DW': np.random.uniform(0, 100),
        
        'Stainless_Steel_p_ALC': np.random.uniform(0, 10),
        'Stainless_Steel_p_CON': np.random.uniform(0, 10),
        'Stainless_Steel_p_DW': np.random.uniform(0, 10),
        'Stainless_Steel_d_ALC': np.random.uniform(0, 10),
        'Stainless_Steel_d_CON': np.random.uniform(0, 10),
        'Stainless_Steel_d_DW': np.random.uniform(0, 10),

        'Stainless_Steel_K_QS_ALC': np.random.uniform(0, 1000),
        'Stainless_Steel_K_QS_CON': np.random.uniform(0, 1000),
        'Stainless_Steel_K_QS_DW': np.random.uniform(0, 1000),

        'Stainless_Steel_A0_ALC': np.random.uniform(0, 100),
        'Stainless_Steel_A0_CON': np.random.uniform(0, 100),
        'Stainless_Steel_A0_DW': np.random.uniform(0, 100),
        
 
    }
    # Different ranges for C0 for KYDEX and Stainless Steel
    for surface in ['KYDEX', 'Stainless_Steel']:
        for cleaning in ['ALC', 'CON', 'DW']:
            if surface == 'KYDEX':
                priors[surface + '_C0_' + cleaning] = np.random.uniform(9000, 15000)
            else:
                priors[surface + '_C0_' + cleaning] = np.random.uniform(5000, 9000)
    
    return priors

def refined_sample_from_prior1():
    priors = {
        # KYDEX ALC
        'KYDEX_r_ALC': np.random.uniform(-2, 2),
        'KYDEX_K_ALC': np.random.uniform(1e4, 3E4),
        'KYDEX_m_ALC': np.random.uniform(5, 10),
        'KYDEX_g_ALC': np.random.uniform(20, 25),
        
        # KYDEX DW
        'KYDEX_r_DW': np.random.uniform(-2, 1),
        'KYDEX_K_DW': np.random.uniform(1E4, 3E4),
        'KYDEX_m_DW': np.random.uniform(7, 10),
        'KYDEX_g_DW': np.random.uniform(17, 24),
        
        # Stainless Steel ALC
        'Stainless_Steel_r_ALC': np.random.uniform(0, 70),
        'Stainless_Steel_K_ALC': np.random.uniform(1E4, 3E4),
        'Stainless_Steel_m_ALC': np.random.uniform(40, 60),
        'Stainless_Steel_g_ALC': np.random.uniform(10, 30),
        
        # Stainless Steel DW
        'Stainless_Steel_r_DW': np.random.uniform(-2, 1),
        'Stainless_Steel_K_DW': np.random.uniform(1E4, 3E4),
        'Stainless_Steel_m_DW': np.random.uniform(2, 5),
        'Stainless_Steel_g_DW': np.random.uniform(15, 20),
    }
    
    # C0 values
    for surface in ['KYDEX', 'Stainless_Steel']:
        if surface == 'KYDEX':
            priors[surface + '_C0'] = np.random.uniform(9000, 25000)
        else:
            priors[surface + '_C0'] = np.random.uniform(5000, 25000)
    
    return priors

def refined_sample_from_prior():
    priors = {
        # KYDEX ALC
        'KYDEX_r_ALC': -0.0262,
        'KYDEX_K_ALC': 19978.4638,
        'KYDEX_m_ALC': 6.7660,
        'KYDEX_g_ALC': 22.4078,
        
        # KYDEX DW
        'KYDEX_r_DW': 0.0905,
        'KYDEX_K_DW': 19997.7876,
        'KYDEX_m_DW': 8.3243,
        'KYDEX_g_DW': 19.5807,
        
        # Stainless Steel ALC
        'Stainless_Steel_r_ALC': 53.2918,
        'Stainless_Steel_K_ALC': 20081.1479,
        'Stainless_Steel_m_ALC': 52.5435,
        'Stainless_Steel_g_ALC': -0.0195,
        
        # Stainless Steel DW
        'Stainless_Steel_r_DW': -0.1978,
        'Stainless_Steel_K_DW': 19999.8625,
        'Stainless_Steel_m_DW': 3.7662,
        'Stainless_Steel_g_DW': 17.8930,
    }
    
    # C0 values
    for surface in ['KYDEX', 'Stainless_Steel']:
        if surface == 'KYDEX':
            priors[surface + '_C0'] = np.random.uniform(8000, 25000)
        else:
            priors[surface + '_C0'] = np.random.uniform(5000, 25000)
    
    return priors



### Distance Functions

def distance_manhattan(simulated_data, observed_data):
    warnings.filterwarnings('ignore', 'You are merging on int and float columns where the float values are not equal to their int representation.')
    warnings.simplefilter(action='ignore', category=UserWarning)
    # Merge simulated and observed data on common keys (Time, Surface, Cleaning)
    merged_data = pd.merge(simulated_data, observed_data, on=['Time', 'Surface', 'Cleaning'])
    # Compute Manhattan differences
    squared_diffs = np.abs(merged_data['mean_sim'] - merged_data['mean'])
    return sum(squared_diffs)


def distance_std(simulated_data, observed_data, epsilon=1e-2):
    merged_data = pd.merge(simulated_data, observed_data.rename(columns={'mean': 'mean_obs', 'std': 'std_obs'}), on=['Time', 'Surface', 'Cleaning'])
    weighted_diffs = np.abs((merged_data['mean_sim'] - merged_data['mean_obs'])) / (merged_data['std_obs'] + epsilon)
    return sum(weighted_diffs)
    



def distance_pcmad(simulated_data, observed_data):
    warnings.filterwarnings('ignore', 'You are merging on int and float columns where the float values are not equal to their int representation.')
    warnings.simplefilter(action='ignore', category=UserWarning)
    # Merge simulated and observed data on common keys (Time, Surface, Cleaning)
    merged_data = pd.merge(simulated_data, observed_data, on=['Time', 'Surface', 'Cleaning'])
    
    # Compute the absolute deviations
    absolute_deviations = np.abs(merged_data['mean_sim'] - merged_data['mean'])
    #print('Abs:',absolute_deviations)
    # Compute the MAD
    mad = median_abs_deviation(absolute_deviations, scale='normal')
    #print('mad:',mad)
    # Identify severe outliers, e.g., those more than 3 MADs from the median
    severe_outliers = absolute_deviations > 3 * mad
    #print('severe:',severe_outliers)
    # Compute the MADO
    mado = median_abs_deviation(absolute_deviations - merged_data['mean'], scale='normal')
    #print('mado:',mado)
    # Compute the PCMAD
    if np.sum(severe_outliers) / len(absolute_deviations) <= 1/3:
        scale = mad
    else:
        scale = mad + mado

    # Compute the scaled absolute differences
    scaled_diffs = absolute_deviations / scale

    return np.sum(scaled_diffs)

#####

## Function that calls the model to fit
def single_simulation(model_func, times, epsilon, observed_data, distance_func):
    params = refined_sample_from_prior()
    simulated_data = simulate_model(model_func, params, times)
    
    if simulated_data['mean_sim'].isnull().any():
        return None

    dist = distance_func(simulated_data, observed_data)
    if dist < epsilon:
        return params, dist  # Return both the parameters and the distance
    return None



def simulate_model(model_func, params, times):
    simulated_data = []
    for surface in ['KYDEX']:#, 'Stainless Steel']:
        for cleaning in ['DW']:#,  'DW']: #'CON',
            
            mean_values = model_func(params, times, surface, cleaning)
            for time, mean in zip(times, mean_values):
                simulated_data.append((time, surface, cleaning, mean))
    return pd.DataFrame(simulated_data, columns=['Time', 'Surface', 'Cleaning', 'mean_sim'])

