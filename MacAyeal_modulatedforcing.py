import matplotlib.pyplot as plt
import numpy as np
import math as m
import matplotlib.gridspec as gridspec
from RK4_func import *

# time steps defined here
# total time = nt * nloop * dt
dt = 0.005 # time step for integration
nt = 8001 # total number of time points
nloop = 100 # number of time points between write out

time = np.linspace(0, nt*nloop*dt, nt) # time array

# time init values
t = 0
t_alpha = 0
t_ramp = 0

x = 1 # initial value for x

# flags initialised here
abateflag = 0
tipflag = 0
incusp = 0
declineflag = 0

switch = 1 # determines if alpha grows (1) or decays (0)

# epsilon == coefficient of the cubic term
epsilon = 0.01 

# solar insolation parameters
qmag = 10
qshift = -6 
T_prec = 20
T_obl = 41
omega_pr = 2*m.pi/T_prec
omega_ob = 2*m.pi/T_obl
A_m = 16

# alpha parameters
alpha_1 = 4.3
alpha_2 = -0.1 
alphareset = -0.05
lambd_a = 0.01 
A_1 = -5.5
A_2 = 4.5

# for plotting purposes, the extra loop is added to make save multiple sims for plotting as subplots on the same fig
plot_x = []
plot_r = []
plot_q = []
plot_alpha = []

# timescale ratios in for the first and second parts of the simulation - Early Pleistocene (EP) and Late Pleistocene (LP)
R_lp_array = [6, 0.9, 25]
R_ep = 25
ramptime = 1 # ramp time (100kyr)


for i in range(0, len(R_mini)):
    R_lp = R_lp_array[i]
    
    # create arrays to store the time series
    R_time = np.zeros(nt)
    alphaarray = np.zeros(nt)
    qarray = np.zeros(nt)
    xx = np.zeros(nt) # time array for plotting the time series of R
    
    t = 0
    t_alpha = 0
    t_ramp = 0
    abateflag = 0
    tipflag = 0
    incusp = 0
    declineflag = 0
    
    for k in range(0, nt):
        # the timescale ratio, R is set here
        if k < 500: 
            R = R_ep
            alpha_1 = 4.3 # 4
            c1 = 1.8 # sinusoidal forcing used with frequency on obliquity timescales
            c2 = 0
        elif k > 501: 
            if t_ramp < ramptime: 
                t_ramp += dt
                c1 = Ramp(1.8, 0, ramptime, t_ramp) # amplitude of sinusoidal obliquity forcing is ramped down to zero
                c2 = Ramp(0, 1, ramptime, t_ramp) # amplitude of precession forcing modulated by eccentricity is ramped up
            else:
                R = R_lp
                alpha_1 = 4.3 # 4
                c1 = 0 # solar forcing on precession timescales modulated by eccentricity is used
                c2 = 1
        for l in range(0,nloop):
            t_alpha += dt
            t += dt
            Q = combined_forcing(A_m, qmag, qshift, omega_ob, omega_pr, t, c1, c2)
            if abateflag == 0:
                alpha = Alphafun(alpha_1, alpha_2, lambd_a, A_1, A_2, R, switch, t_alpha)
            else:
                alpha = alphareset - 0.01
            if alpha>0:
                cusp = CuspLocus(alpha, epsilon)
            
            x = RK4(x, epsilon, t, alpha, Q, dt) # integrate
 
            # if the system in onset phase and in the bistable region, set incusp flag
            if declineflag==0 and alpha>0 and Q > cusp*-1:
                 incusp=1
            
            # if the system has passed over a saddle node bifurcation, set reverse tip flags, reset other flags
            if tipflag==1:
                incusp=0
                tipflag=0
            
            # if the incusp flag is 1, the system is now outside of the cusp locus, and a bifurcation is taking place
            # alpha switches from growth to decay
            if incusp==1 and tipflag==0 and alpha>0 and Q > cusp:
                tipflag = 1
                declineflag = 1
                switch = 0
                A_2 = abs(alpha_2-alpha)
                t_alpha = 0
                            
            # if alpha is decaying then the system should not switch to growth until alpha is less than the reset value
            # if Q is less than 0 then the alpha onset is aborted until abate flag is == 0
            if declineflag==1 and alpha<alphareset: 
                if Q > 0:
                    abateflag = 1
                    # solar insolation is too high for onset
                    # alpha remains negative
                else:
                    # alpha onset can occur
                    # switch to alpha growth
                    abateflag = 0
                    declineflag = 0
                    switch = 1
                    A_1 = -1*abs(alpha_1-alpha)
                    t_alpha = 0
        
        x_saved = x
        xx[k] = x
        alphaarray[k] = alpha
        qarray[k] = Q
        R_time[k] = R
    
    plot_x.append(xx*-1)
    plot_r.append(R_time)
    plot_alpha.append(alphaarray)
    plot_q.append(qarray)
    
# End of loop


# ALL PLOTTING BELOW

# plot the cusp locus and the trajectory in parameter space
locus_pos = np.zeros(100)
locus_neg = np.zeros(100)
alphalocus = np.linspace(0,3,100)

for i in range(0, 100):
    locus_pos[i] = CuspLocus(alphalocus[i], epsilon)
    locus_neg[i] = -1*CuspLocus(alphalocus[i], epsilon)

fig1 = plt.figure()
plt.plot(qarray, alphaarray, 'black')
plt.plot(locus_pos, alphalocus, 'r')
plt.plot(locus_neg, alphalocus, 'r')
plt.xlabel('Solar insolation')
plt.ylabel('Alpha')


x = np.linspace(-40, 40, 400)
y = np.linspace(-5, 5, 400)
x, y = np.meshgrid(x, y)
z = (epsilon*(x**3) - x*y)

y_line = np.linspace(1, 10, 400)
z_line = np.sqrt((2*y_line**3)/(27*epsilon))

fig3 = plt.figure()
ax = fig3.add_subplot(111, projection='3d')
ax.plot_surface(z, y, x, cmap='Greys', alpha=0.4)
ax.set_xlim([-20,20])
ax.set_ylim([-5,5])
ax.set_zlim([-20,20])
ax.auto_scale_xyz([-20,20], [0,5], [-20,20])
ax.set_xlabel("Q")
ax.set_ylabel("Alpha")
ax.set_zlabel("X")

ax.plot(z_line, y_line, 'black')
ax.plot(-1*z_line, y_line, 'black')
ax.plot(plot_q[1][2000:-1], plot_alpha[1][2000:-1], plot_x[1][2000:-1]*-1, 'red')



fig = plt.figure(figsize=(20, 6))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.1, wspace= 0.04, height_ratios=[2, 2, 1])

with open('Documents/Time_ka_Benthic_d18O.txt', 'r') as file:
# Read the lines into an array
    data = file.readlines()

# Strip newline characters from each line
data = [line.replace('\t', ',') for line in data]
data = [line.strip('\n') for line in data]

# Function to safely convert strings to floats
def safe_convert(s):
    try:
        return float(s)
    except ValueError:
        return 'Hi'
        
# Convert each string to a list of floats, handling potential issues
float_array = []
switch1 = 0
for s in data:
# Split the string by commas
    split_values = s.split(',')
    # Convert each value to float if possible, otherwise keep None
    #float_row = [safe_convert(value) for value in split_values]
    for value in split_values:
        val = safe_convert(value)
        if type(val) == str:
            switch1=1
    float_row = [safe_convert(value) for value in split_values]
    if switch1==0:
        float_array.append(float_row)
    switch1=0
    
reshaped_array = np.array(float_array).T
data = reshaped_array
y_top = np.full(100, 2.7)
y_bottom = np.full(100, 5.3)
x_line = np.linspace(800, 1200, 100)


ax = fig.add_subplot(gs[0, 0:3])
color = 'black'
width = 0.9
ax.plot(data[0], data[1], color=color, linewidth=width)
ax.invert_xaxis()
ax.xaxis.tick_top()
ax.tick_params(axis="x", direction="in")

ax.fill_between(x_line, y_top, y_bottom, alpha=0.3, color='m')
ax.set_ylabel(r'$\delta$ O$^{18}$ ($\text{‰}$)')
ax.set_xlabel('Time (kyr)', fontsize=15)
ax.xaxis.set_label_position('top')
ax.set_xticks(np.linspace(2000, 0, 11))
ax.text(x=1000, y=2.8, s="MPT", ha='center', va='bottom')
ax.set_xlim(2000, 0)
ax.set_ylim(2.7, 5.3)
gs.update(hspace=0.1)


ax1 = fig.add_subplot(gs[1, 1:3])
color = 'black'
width = 0.9
ax1.plot(time, plot_x[2], color=color, linewidth=width)
ax1.set_ylabel('State Variable')
ax1.set_xlim(0, 1600)
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax1.yaxis.set_label_position('right')
ax1.yaxis.tick_right()

color = 'dodgerblue'
ax2 = fig.add_subplot(gs[2, 1:3])
ax2.plot(time, plot_alpha[2], 'dodgerblue', linewidth=width)
ax2.set_ylabel('Alpha', labelpad=5)
ax2.tick_params(axis='y')
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()
ax2.set_xlim(0, 1600)
ax2.set_ylim(-0.5, 3.5)
ax2.set_xlabel('Time (kyr)')

fig.align_ylabels([ax1, ax2])

ax = fig.add_subplot(gs[1:3, 0])
color = 'black'
ax.plot(plot_q[2][2000:-1], plot_alpha[2][2000:-1], color=color, linewidth=width)
ax.set_ylim(-0.5, 3)
ax.plot(locus_pos, alphalocus, color='red', linewidth=width, linestyle='--')
ax.plot(locus_neg, alphalocus, color='red', linewidth=width, linestyle='--')
ax.tick_params(axis='y', labelcolor=color)
ax.tick_params(axis="y", direction="in")
ax.tick_params(axis='x', labelcolor=color)
ax.set_ylabel('Alpha', labelpad=0)#rotation=270, labelpad=10)
ax.set_xlabel('Solar Insolation', labelpad=0)

plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

plt.show()



# Lisiecki, L. E., and Raymo,  M. E. (2005) A Pliocene-Pleistocene stack of 57 globally distributed benthic d18O records, Paleoceanography, 20, PA1003, doi:10.1029/2004PA001071
