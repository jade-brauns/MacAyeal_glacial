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
epsilon = 0.02 

# solar insolation parameters
q1 = 6
q0 = 4
T_q = 41 # period of q (ka)
omega_q = 2*m.pi/T_q

# alpha parameters
alpha_1 = 4 
alpha_2 = -0.5 
alphareset = -0.4 
lambd_a = 0.04 # the decay constant for alpha onset, while alpha decline had decay constant R*lamb
A_1 = -5.5
A_2 = 4.5

# initialize arrays for storing simulated data
plot_x = []
plot_r = []
plot_q = []
plot_alpha  = []

# timescale ratios in for the first and second parts of the simulation - Early Pleistocene (EP) and Late Pleistocene (LP)
R_ep = 4
R_lp_array = [3, 0.8, 0.6]




# for plotting purposes, the extra loop is added to make save multiple sims for plotting as subplots on the same fig
for i in range(0, len(R_lp_array)):
    R_lp = R_lp_array[i]
    ramptime = 1 # ramping time is 100 kyr
    # create arrays to store the time series
    R_time = np.zeros(nt) # time array for plotting the time series of R
    alphaarray = np.zeros(nt)
    qarray = np.zeros(nt)
    xx = np.zeros(nt)
    
    t = 0
    t_alpha = 0
    t_ramp = 0
    abateflag = 0
    tipflag = 0
    incusp = 0
    declineflag = 0
    
    for k in range(0, nt):
        # the timescale ratio, R is set here
        if k < 401: 
            R = R_ep
        elif k > 400:
            if t_ramp < ramptime: # R is linearly ramped down
                t_ramp += dt
                R = Ramp(R_ep, R_lp, ramptime, t_ramp)
            else:
                R = R_lp
                
        for l in range(0,nloop):
            t_alpha += dt
            t += dt
            Q = Solarforcing(q1, q0, omega_q, t)
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
                    # solar insolation is too high for alpha onset
                    # alpha remains negative
                else:
                    # alpha onset can occur
                    # switch to alpha growth
                    abateflag = 0
                    declineflag = 0
                    switch = 1
                    A_1 = -1*abs(alpha_1-alpha)
                    t_alpha = 0
        
        xx[k] = x
        alphaarray[k] = alpha
        qarray[k] = Q
        R_time[k] = R
    
    plot_x.append(xx*-1) # glacial volume is equal to -x following the convention of MacAyeal (1979)
    plot_r.append(R_time)
    plot_alpha.append(alphaarray)
    plot_q.append(qarray)
    
# End of loop


# ALL PLOTTING BELOW

# plot the cusp locus and the trajectory in parameter space
locus_pos = np.zeros(100)
locus_neg = np.zeros(100)
alphalocus = np.linspace(0,2.2,100)

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

fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
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


fig = plt.figure(figsize=(20, 14))
gs = gridspec.GridSpec(6, 3, figure=fig, hspace=0.1, wspace= 0.15, height_ratios=[2, 1, 2, 1, 2, 1])
gs.update(hspace=0.1)


count=0
count1=0
# Plotting the remaining three subplots in separate rows
for i in range(3):

    ax1 = fig.add_subplot(gs[count, 1:3])
    color = 'black'
    width = 0.9
    ax1.plot(time, plot_x[i], color=color, linewidth=width)
    if i==1:
        ax1.set_ylabel('State Variable')
    ax1.set_xlim(0, 1500)
    ax1.set_ylim(-12, 10)
    ax1.tick_params(axis='y', labelcolor=color)
    ax_twin = ax1.twinx()
    color = 'tab:red'
    ax_twin.plot(time, plot_r[i], color=color, linewidth=1)
    if i==1:
        ax_twin.set_ylabel('Timescale Separation', color=color)
    ax_twin.yaxis.set_label_position('right')
    ax_twin.tick_params(axis='y', labelcolor=color)

    # Hide ticks for the combined plots
    if count<6:
        ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax_twin.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    if i<4:
        count+=1
        color = 'royalblue'
        ax2 = fig.add_subplot(gs[count, 1:3])
        ax2.plot(time, plot_alpha[count1], 'royalblue', linewidth=width)
        if i==1:
            ax2.set_ylabel('Alpha', color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_xlim(0, 1500)
        ax2.set_ylim(-0.5, 2.3)
        if i==2:
            ax2.set_xlabel('Time (Kyr)')
            fig.align_ylabels([ax1, ax2])
        count1+=1
    if count<5:
        ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax_twin.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    if i==1:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax_twin.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    count+=1

    

# Plotting the three subplots in the final column (state space)
for i in range(0,3):
    ax = fig.add_subplot(gs[(2*i):(2*i+2), 0])
    color = 'black'
    ax.plot(plot_q[i][2000:-1], plot_alpha[i][2000:-1], color=color, linewidth=width)
    ax.set_ylim(-0.5, 2.2)
    ax.plot(locus_pos, alphalocus, color='red', linewidth=width, linestyle='--')
    ax.plot(locus_neg, alphalocus, color='red', linewidth=width, linestyle='--')
    ax.tick_params(axis='y', labelcolor=color)
    ax.tick_params(axis="y", direction="in")
    ax.tick_params(axis='x', labelcolor=color)
    ax.set_xlim(-5, 11)
    if i==0:
        newax = fig.add_axes([0.125,0.82,0.05,0.06], zorder=5)
        newax.plot(plot_q[i][100:500], plot_alpha[0][100:500], color=color, linewidth=width)
        newax.plot(locus_pos, alphalocus, color='red', linewidth=width, linestyle='--')
        newax.plot(locus_neg, alphalocus, color='red', linewidth=width, linestyle='--')
        newax.set_xlim(-5, 11)
        newax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        newax.set_yticks([])        
    if i<2:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        if i==1:
            ax.yaxis.set_label_position('left')
            ax.set_ylabel('Alpha', labelpad=0)#rotation=270, labelpad=10)
    else:
        ax.set_xlabel('Q')


plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
plt.show()

