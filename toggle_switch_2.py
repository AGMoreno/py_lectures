# Create toggle switch model
from scipy.optimize import fsolve
from scipy.integrate import odeint, quad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# ODEs 
def f_ode(x,t, *args):
    K,b_t,b_l,d_t,d_l,del_t,del_l,beta_t,beta_l = args
    m_t = x[0]
    m_l = x[1]
    p_t = x[2]
    p_l = x[3]
    y = [0,0,0,0]
    y[0] = K * b_t/(b_t + p_l) - d_t * m_t
    y[1] = K * b_l/(b_l + p_t) - d_l * m_l
    y[2] = beta_t * m_t - del_t*p_t
    y[3] = beta_l * m_l - del_l * p_l
    return y

x0 = [10,5,50,0]
params_nom = [100, 100, 10, 5, 5, 0.01, 0.01, 0.01, 0.01]
# params_nom = [0.01, 0.01, 0.01, 1, 1, 0.1, 0.1, 1, 1]

timepoints = np.linspace(0,0.5,100)
params = tuple(params_nom)
# Solve ODEs, numerically
y = odeint(f_ode, x0, timepoints, args = params)

# Calculate Lipscitz constant
K,b_t,b_l,d_t,d_l,del_t,del_l,beta_t,beta_l = params_nom
L = np.max([del_l + K/b_t, beta_l + d_t, beta_t + d_l, del_t + K/b_l])
print('The Lipschtiz constant is {0}'.format(L))
# Check if the L holds globally

# Time varying bound
# norm(xt - zt) <= norm(x0 - z0) exp(L*(t-t0)) + mu * (t-t0) * exp(L*(t-t0))
xtzt = np.zeros((len(timepoints),4))
t0 = 0
# Now vary parameters to find y_zt, then plot y-y_zt
params_new = [i for i in params_nom]
# K perturbed by Delta = 10 (K = K + Delta)
Delta = 0
params_new[0] += Delta
# According to the dynamics, we can show that this corresponds to the following mu
mu = 2*Delta
z0 = [i+5 for i in x0]
for i in range(len(timepoints)):
    t = timepoints[i]
    for s in range(4):
        diff_ic = x0[s] - z0[s]
        xtzt[i,s] = np.linalg.norm(diff_ic) * np.exp(L*(t-t0)) + mu * (t-t0) * np.exp(L*(t-t0))

y_zt = odeint(f_ode, x0, timepoints, args = tuple(params_new))

y_err = np.zeros((len(timepoints),4))
for s in range(4):
    for t in range(len(timepoints)):
        y_err[t][s] = y[t][s] - y_zt[t][s]

for i in range(4):
        plt.subplot(2,2,i+1)
        # Plot bound
        plt.plot(timepoints, xtzt[:,i],'--',linewidth = 1.5, label = 'bound')
        # Plot solutions from numerical simulation
        plt.plot(timepoints, y_err[:,i], label = 'x(t) - z(t)')
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('State x_' + str(i+1))
plt.show()