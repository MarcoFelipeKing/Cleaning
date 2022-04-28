# README

This is the cleaning git readme. Each sub-folder contains a README of its own. But for brevity, they are pasted in here.


To do:
Find mean and standard deviation of stephanie's data based on cleaning method.
Fit ODEs to Stephanie's data to find parameters for the PAMCODE

-----

# StephanieData


Bogusz 2013 contains the detergent based cleaning data over 4 months (once per month)
120 surfaces (30 rooms, 4 surfaces per room: left and right bed rails, beside locker and overbed table). They use detergent Tuffie wipes

Stewart 2014 is the 2013 experiment repeated with electrolysed water.


allCleanDataStephanie is the amalgamation of the results for both studies. It's cbased on her categories: no growth 0/cm2, scanty growth >2.5cfu/cm2, low growth >12/cm2, medium growth >26/cm2 and high growth >40/cm2. These numbers might not be right, but they're an approximation from memory (on 12th of August 2020). Data in the file is faithfully transcribed though.

##########
# LabData

Inoculated (3*3cm) surfaces with S. aureus. Allowed to dry for 60 min.
Cleaned with either Disinfectant wipe (Tuffie), detergent wipe (Tuffie), distilled water wipe (Tuffie) or no cleaning=control
Swabbed every hour for 8h then at 24h.
5 concurrent replicates.


##########
# CODE

######## ABC ###########
# Update 12th August

ODE parameters for the Alcohol curve are in parameter_sample_Alcohol_1K.csv. These were obtained using the pyABC_Cleaning.py code.
So technically we can load these values directly into the PAM_Cleaning_Curves.R code instead of the random values we have currently.

Cleaning_ABC.R is an ABC rejection algorithm in R. There is also an SMC algorithm undereath it but I can't get it to work. I have a python code that works fine though. (pyABC_Cleaning.py)

To run Cleaning_ABC.R all you need to to is change the experimental data hard-coded into it. it will spit out a dataframe with 1000 best values for the parameters.
initial_contamtination=c(Contamination=59) #This is the initial value immediately after cleaning (time=0)
experimental_data= c(Contamination@time1,Contamination@time2,...,Contamination@time24)
s=standard deviation

To run the pyhon code, I recommend installing anaconda and using spyder. I've had no end of hassle with stand-along python installations.



######## PAMCODE ###########

This is based on the PPE code but includes the cleaning curve ODES. as yet it uses random values for the parameters, these need estimating first using either the Cleaning_ABC.R code or the pyABC_Cleaning.py SMC code. Should just run and produce boxplots of exposure at each time-point. Bit rough around the edges (assumes only 1 hand is used and homogeneous spread around the room)
