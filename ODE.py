import scipy
from scipy import integrate
from pylab import *
from operator import itemgetter

v = scipy.zeros((4),'d')
def rate_eqs(S):
    mRNA=S[0];  gfp=S[1]; 
    v[0] = k1 * P
    v[1] = k4 * gfp
    v[2] = k3 * mRNA
    v[3] = k2 * mRNA
    return v

def diff_eqs(S,t):
    Sdot = scipy.zeros((2),'d')
    v = rate_eqs(S)
    Sdot[0] =  + v[0] - v[2]
    Sdot[1] =  - v[1] + v[3]
    return Sdot

mRNA=0.0;  gfp=0.0;  P=18.0; 
k1=1.5;  k4=0.265;  k3=0.65;  k2=2.5; 

t_start=0.0; t_end=100.0; t_inc=0.1
t_range = scipy.arange(t_start, t_end+t_inc, t_inc)
t_course = scipy.integrate.odeint(diff_eqs, [mRNA,gfp], t_range)

plot(t_range, map(itemgetter(0),t_course))
plot(t_range, map(itemgetter(1),t_course))
show()

