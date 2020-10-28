import numpy as np
import random as r
import math as m
import matplotlib.pyplot as plt

mRNA=0.0;  gfp=0.0;  P=18.0; 
k1=1.5;  k4=0.265;  k3=0.65;  k2=2.5; 
t=0.0
tmax=500.0

trange = []
mRNArange = []
gfprange = []
Prange = []

while t < tmax:
    a = [ k1 * P , k4 * gfp , k3 * mRNA , k2 * mRNA ]
    a0 = sum(a)
    r1 = r.random()
    tau = -m.log(r1)/a0
    t = t+tau
    r2 = r.random()
    acumsum = np.cumsum(a)/a0
    chosen_reaction = min([i for i in range(len(a)) if acumsum[i] >= r2])
    if chosen_reaction ==  0 :
         mRNA+=1;
    if chosen_reaction ==  1 :
         gfp-=1;
    if chosen_reaction ==  2 :
         mRNA-=1;
    if chosen_reaction ==  3 :
         gfp+=1;
    trange.append(t)
    mRNArange.append(mRNA)
    gfprange.append(gfp)
    Prange.append(P)

plt.plot(trange,mRNArange)
plt.plot(trange,gfprange)
plt.plot(trange,Prange)
plt.show()

