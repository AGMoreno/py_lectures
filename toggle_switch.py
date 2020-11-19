# Create toggle switch model
# https://github.com/ayush9pandey/toggle_switch/blob/master/Project%20presentation.pdf
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
    y[0] = K * b_t**2/(b_t**2 + p_l**2) - d_t * m_t
    y[1] = K * b_l**2/(b_l**2 + p_t**2) - d_l * m_l
    y[2] = beta_t * m_t - del_t*p_t
    y[3] = beta_l * m_l - del_l * p_l
    return y

x0 = [0,0,0,0]
params_nom = [100, 10, 100, 5, 5, 0.01, 0.01, 0.01, 0.01]
params = tuple(params_nom)

timepoints = np.linspace(0,500,100)
# Solve ODEs, numerically for different ic
y = odeint(f_ode, x0, timepoints, args = params)
plt.plot(timepoints, y[:,0],linewidth = 1.5,label='mRNA TetR')
plt.plot(timepoints, y[:,1],linewidth = 1.5, label='mRNA LacI')
plt.plot(timepoints, y[:,2],linewidth = 1.5, label = 'protein TetR')
plt.plot(timepoints, y[:,3],linewidth = 1.5, label = 'protein LacI')
plt.xlabel('Time')
plt.ylabel('Amount')
plt.legend()
plt.show()
