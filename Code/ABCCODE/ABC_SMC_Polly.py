import numpy as np
from scipy import integrate

#import sys
#index = sys.argv[1]

spore_data_mean = np.array((59,19,5.4,5.4,2.4,8.6))#pd.readcsv('../../LabData/cleaningExptBeth_summary.csv')
spore_data_sd = np.array((26,2.3,4.7,4.3,4.3))#pd.readcsv('Spore_data_sd.txt')
#bac_data_mean = np.loadtxt('Bacterial_data_mean.txt')
#bac_data_sd = np.loadtxt('Bacterial_data_sd.txt')

#ICs = #np.loadtxt('ICs.txt')
s00 = np.array((0))#ICs[0]
s10 = np.array((0))#ICs[1]

#Define the ordinary differential equation models 
def dN_dt(N,t,par):
   C,r,g,m = par.T
   return np.array([
          r*N[0]*(1-N[0]/C)-m*np.exp(-g*t),
         ])   
   
#Set the timecourse for the integration (hours)
tt = np.linspace(0.0, 24, 501)

#Function to generate a first value of epsilon to use in the ABC
def first_iteration(N,ep):
   print('Running iteration: 0')
   epsilon = ep
   #Empty arrays for the accepted parameters of each hypothesis
   accepted_params = np.empty((0,6))
   #Store the distance measures from the ABC
   results = np.empty((0))
   number = 0 #Counter for the population 
   truns = 0 #Count total number of runs to get acceptance percentage
   while number < N:
      truns+=1
      #Draw the parameters from the prior distributions if it is the first iteration
      C = 10**np.random.normal(0,100)
      r = 10**np.random.uniform(0,1)
      g = 10**np.random.uniform(0,1)#np.random.normal(np.log10(0.8),0.5)
      m = 10**np.random.uniform(0,1)
           
      eps1 = np.random.uniform(0,1)     
      eps2 = np.random.uniform(0,1)     
      
      s0a = s00*eps1
      s0b = s00*(1-eps1)
      s1a = s10*eps2
      s1b = s10*(1-eps2)
  
      parset = np.array((C,r,g,m)) #Parameter set 
      parset2 = np.array((C,r,g,m,eps1,eps2))
      R0 = np.array((59.4)) #Initial concentrations
      R = integrate.odeint(dN_dt, R0, tt, args=(parset,)) #Integrate the mathematical model
      n1 = R.T #Model outputs for each variable
      spore = n1 + s0b #Total model CFU
      #bac = n2 + n3

      #Take the model points of the model simulations that correlate to the data time points
      time_points = np.array([0,1,2,4,8,24])
      spore_output = np.log10(spore[time_points])
      #bac_output = np.log10(bac[time_points])
      
      spore_final = spore_output + np.random.normal(0,spore_data_sd)
      #bac_final = bac_output + np.random.normal(0,bac_data_sd)

      dists_spore = np.sum(pow((spore_final - spore_data_mean),2))
      #dists_bac = np.sum(pow((bac_final - bac_data_mean),2))
      distance = np.sqrt(dists_spore)#np.sqrt(dists_spore + dists_bac)

      if distance < epsilon:
         number+=1
         results = np.hstack((results,distance))
         accepted_params = np.vstack((accepted_params, parset2))
            
   #Weights
   weights = np.empty((0,1))
   for i in range(len(accepted_params)):
      weights = np.vstack((weights,1/len(accepted_params)))
      
   print('Acceptance rate for iteration 0: ' + str(N*100/truns))
   print('Epsilon = ' + str(epsilon))
   print('Total runs = ' + str(truns))
   return [np.hstack((np.reshape(results,(len(accepted_params),1)),accepted_params,weights)),truns]

def other_iterations(N,it):
   print('Running iteration: ' + str(it+1))
   
   epsilon = np.median(ABC_runs[it][:,0])
   #epsilon = epsilons[it]
   p_list = [i for i in range(N)]

   #Upper and lower bounds for the uniform distributions of the priors
   lower_bounds = [0, 0, 0, 0, 0, 0]
   upper_bounds = [100, 10, 10, 10, 10, 10] 
   
   #Uniform areas to sample within in order to perturb the parameters
   ranges = []
   for i in range(6):
      if i in [4,5]:
         r1 = np.max(ABC_runs[it][:,i+1]) - np.min(ABC_runs[it][:,i+1])
      else:
         r1 = np.max(np.log10(ABC_runs[it][:,i+1])) - np.min(np.log10(ABC_runs[it][:,i+1]))
      ranges.append(r1)
   ranges_arr = np.asarray(ranges)
   
   sigma = 0.2*ranges_arr

   #Empty arrays for the prior samples
   priors = np.empty((0,6))
   #Empty arrays for the accepted parameters of each hypothesis
   accepted_params = np.empty((0,6))
   #Store the distance measures from the ABC
   results = np.empty((0))
   #Store the weights for each hypothesis
   weights = np.empty((0))

   number = 0
   truns = 0
   while number < N:
      #print(number)
      truns+=1
      check = 0
      while check < 1:
         choice = np.random.choice(p_list,1,p=ABC_runs[it][:,7])
         prior_sample = ABC_runs[it][:,range(1,7)][choice]
         parameters = []
         for i in range(6):
            if i in [4,5]:
               lower = prior_sample[0,i]-sigma[i]
               upper = prior_sample[0,i]+sigma[i]
               parameter = np.random.uniform(lower, upper)
               parameters.append(parameter)
            else:
               lower = np.log10(prior_sample[0,i])-sigma[i]
               upper = np.log10(prior_sample[0,i])+sigma[i]
               parameter = np.random.uniform(lower, upper)
               parameters.append(pow(10,parameter))

         check_out = 0
         for ik in range(6):
            if ik in [4,5]:
               if parameters[ik] < lower_bounds[ik] or parameters[ik] > upper_bounds[ik]:
                  check_out = 1
            else:
               if parameters[ik] < pow(10,lower_bounds[ik]) or parameters[ik] > pow(10,upper_bounds[ik]):
                  check_out = 1
         if check_out == 0:
            check+=1

      C = parameters[0]
      r = parameters[1]
      g = parameters[2]
      m = parameters[3]
      eps1 = parameters[4]
      eps2 = parameters[5]

      s0a = s00*eps1
      s0b = s00*(1-eps1)
      s1a = s10*eps2
      s1b = s10*(1-eps2)

      parset = np.array((C,r,g,m)) #Parameter set 
      parset2 = np.array((C,r,g,m,eps1,eps2))
      R0 = np.array(59.4) #Initial concentrations
      R = integrate.odeint(dN_dt, R0, tt, args=(parset,)) #Integrate the mathematical model
      n1 = R.T #Model outputs for each variable
      spore = n1 #+ s0b #Total model CFU
      #bac = n2 + n3

      #Take the model points of the model simulations that correlate to the data time points
      time_points = np.array((0,1,2,4,8,24))
      spore_output = np.log10(spore[time_points])
      #bac_output = np.log10(bac[time_points])

      spore_final = spore_output + np.random.normal(0,spore_data_sd)
      #bac_final = bac_output + np.random.normal(0,bac_data_sd)

      dists_spore = np.sum(pow((spore_final - spore_data_mean),2))
      #dists_bac = np.sum(pow((bac_final - bac_data_mean),2))
      distance = np.sqrt(dists_spore)

      if distance < epsilon:
         number+=1         
         denom_arr = []
         for j in range(N):
            weight = ABC_runs[it][j,7]
            params_row = ABC_runs[it][j,1:7]
            boxs_up = []
            boxs_low = []
            for i in range(6):
               if i in [4,5]:
                  boxs_up.append(params_row[i] + sigma[i])
                  boxs_low.append(params_row[i] - sigma[i])
               else:
                  boxs_up.append(np.log10(params_row[i]) + sigma[i])
                  boxs_low.append(np.log10(params_row[i]) - sigma[i])
            outside = 0
            for i in range(6):
               if i in [4,5]:
                  if parameters[i] < boxs_low[i] or parameters[i] > boxs_up[i]:
                     outside = 1                                    
               else:
                  if np.log10(parameters[i]) < boxs_low[i] or np.log10(parameters[i]) > boxs_up[i]:
                     outside = 1                  
            if outside == 1:
               denom_arr.append(0)
            else:
               denom_arr.append(weight*np.prod(1/(2*sigma)))

         weight_param = 1/np.sum(denom_arr)
         
         weights = np.hstack((weights,weight_param))
         results = np.hstack((results, distance))
         accepted_params = np.vstack((accepted_params, parset2))
         priors = np.vstack((priors, prior_sample))

   weights_2 = weights/np.sum(weights)
   weights_3 = np.reshape(weights_2, (len(weights_2),1))
   
   print('Acceptance rate for iteration ' + str(it+1) + ': ' + str(N*100/truns))
   print('Epsilon = ' + str(epsilon))
   print('Total runs = ' + str(truns))
   return [np.hstack((np.reshape(results,(len(accepted_params),1)),accepted_params,weights_3)),truns]

#Sample size for the ABC
N = 500
first = first_iteration(N,100)
ABC_runs = []
ABC_runs.append(first[0])
np.savetxt('ABC_set_' + str(1) + '_run_0.txt', ABC_runs[0])
#Run the sucessive iterations of the ABC
for itt in range(30):   
   run = other_iterations(N,itt)
   ABC_runs.append(run[0])
   np.savetxt('ABC_set_' + str(1) + '_run_' + str(itt+1) + '.txt', ABC_runs[itt+1])

