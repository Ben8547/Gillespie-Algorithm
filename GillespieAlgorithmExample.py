import numpy as np
import matplotlib.pyplot as plt
import random as rand
from scipy.special import comb

#maxwell-boltzman distrobution

# Actual attempt for DRG project:

#Reaction 29: X + Y -> 2Y (c_1)
#            2Y -> Z (c_2)

# c * the number of particles of reactant 1 * partices of reactant 2 * dt tells us the probabilty that the reaction occurs in t -> t+dt

c_1X = 5 #(We assume that X is in excess so that this value does not change)

c_2 = 0.00125

Y_int_1 = 40

Y_int_2 = 12000

prob_1 = lambda Y: c_1X*Y

prob_2 = lambda Y: c_2*comb(Y,2)

def reaction_prob_density_func(tau,mu,state:list):
    # We assume [X] to be constant to ombit it from "state"; Y should be first entry, and Z should be the second (though Z not needed either).
    if mu == 1:
        a_mu = prob_1(state[0])
        a_0 = a_mu + prob_2(state[0])
    elif mu == 2:
        a_mu = prob_2(state[0])
        a_0 = a_mu + prob_1(state[0])
    else:
        raise Exception("Invalid Reaction channel")

    return a_mu * np.exp(-a_0*tau)

#Step:

def step(r_1,r_2,state):

    #the state variable should be [Y,Z,t]
    
    a_mu_1 = prob_1(state[0])

    a_mu_2 = prob_2(state[0])

    a_0 = prob_2(state[0]) + prob_1(state[0])

    tau = (1/a_0) * np.log(1/r_1)

    mu=[]
    for i in [1,2]:
        if i == 1:
            if 0 < r_2*a_0 and r_2*a_0 <= a_mu_1:
                mu.append(i)
        if i == 2:
            if a_mu_1 < r_2*a_0 and r_2*a_0 <= np.sum([a_mu_1,a_mu_2]):
                mu.append(i)
    #print(mu)

    if mu == []:
        #no reaction occurs - state does not change except in time
        mu=[] # do nothing
    if 1 in mu:
        state[0] += 1
    if 2 in mu:
        state[0] = state[0] - 2
        state[1] += 1
    
    state[2] += tau

    return state

#Plot

init_state_1 = [Y_int_1,0,0]
init_state_2 = [Y_int_2,0,0]

reaction_trial_1 = [Y_int_1]

reaction_trial_2 = [Y_int_2]

time_1 = [0]

time_2 = [0]

while time_1[-1] < 5:
    r_11 = rand.random()
    r_12 = rand.random()
    reaction_trial_1.append((init_state_1:=step(r_11,r_12,init_state_1))[0])
    time_1.append(init_state_1[2])

while time_2[-1] < 5:
    r_21 = rand.random()
    r_22 = rand.random()
    reaction_trial_2.append((init_state_2:=step(r_21,r_22,init_state_2))[0])
    time_2.append(init_state_2[2])

fig, ax = plt.subplots()

ax.plot(time_1, reaction_trial_1, ls='-',marker='none',label = 'Initial Y = 40')

ax.plot(time_2, reaction_trial_2,ls='-',marker='none',label='Initial Y = 12000')

ax.axhline(y=4000,label = "Steady State Equilibrium, Y=4000",color='cyan',ls=':')

ax.set_ylabel("Number of Y molecules")

ax.set_xlabel("Time")

ax.set_title("Reaction of X + Y -> 2Y -> Z with X in large excess")

#print(time_1)
#print(time_2)

ax.legend()

plt.show()
