#script to fit parameters of ODE model to data using ABC
#using DifferentialEquations
using DifferentialEquations
using Distributions

#function for ODE for calculating the number of bacteria on a surface at time t
# r*(1-Bacteria/C)*Bacteria-m*math.exp(-g*t)*Bacteria

function f(du,u,params,t)
    r, C, m, g = params
    du[1] = r*(1-u[1]/C)*u[1]-m*exp(-g*t)*u[1]
end

# function to calculate the distance between the model and the data
function distance(data, model)
    return sum(abs.(data .- model))
end

# function to calculate the model and return bacteria t time points 0, 1, 2, 4, 12 and 24 hours
function model(params)
    u0 = [221.6]
    tspan = (0.0,24.0)
    prob = ODEProblem(f,u0,tspan,params)
    sol = solve(prob,Tsit5(),saveat=0:24)
    #return only the bacteria at the time points 0, 1, 2, 4, 12 and 24 hours
    return sol[1,1:6]
end

# experimental data for bacteria at time points 0, 1, 2, 4, 12 and 24 hours
data_mean =  [221.6,94.3,56.25,1.75,1.6,8.5]
data_sd = [76.4,86.9,89.4,0.5,2.3,4.04]

#set up the prior distribution
#r ~ Uniform(0,1)
#C ~ Uniform(0,1000)
#m ~ Uniform(0,1)
#g ~ Uniform(0,1)
prior = [Uniform(0,1), Uniform(0,1000), Uniform(0,1), Uniform(0,1)]


# for 10000 simulations choose a value for each variable from their prior
# calculate the model
# calculate the distance between the model and the data
# if the distance is less than 10 then save the parameters
# if the distance is greater than 10 then discard the parameters
# repeat until 1000 parameters have been saved

#create and empty array to store the parameters
parameters = []

#set up the ABC algorithm
while true
    #choose a value for each variable from their prior
    r = rand(prior[1])
    C = rand(prior[2])
    m = rand(prior[3])
    g = rand(prior[4])
    params = [r, C, m, g]
    #calculate the model
    model_data = model(params)
    #calculate the distance between the model and the data
    d = distance(data_mean, model_data)
    #if the distance is less than 10 then save the parameters and the distance to an array parameters
    if d < 25
        push!(parameters, [params, d])
    end
    #if the distance is greater than 100 then discard the parameters
    #repeat until 1000 parameters have been saved
    if length(parameters) == 1000
        break
    end
end

#sort the parametrs by distance
parameters = sort(parameters, by = x -> x[2])

using DelimitedFiles
#save the parameters and distances to a file csv file  ABC_parameters
writedlm("ABC_parameters.csv", parameters)

#plot the best parameter combination against the data with standard deviation data_sd
using Plots
plot([0,1,2,4,12,24], model(parameters[1][1]), label = "Model", xlabel = "Time (hours)", ylabel = "Bacteria (CFU)", title = "ABC model fit to data", legend = :topleft)    
#add errorbars to show the standard deviation of the data
scatter!([0,1,2,4,12,24], data_mean, yerr = data_sd, label = "Data", legend = :topright)


