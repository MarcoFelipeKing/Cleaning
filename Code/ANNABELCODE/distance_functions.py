# distance functions

## Distance Functions

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
