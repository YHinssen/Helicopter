from math import *
import numpy as np
import matplotlib.pyplot as plt

# ------ Constants ---------#
k = 1.15
W = 93440.25 #N
rho = 1.225
R = 7.315 #m
N_blades = 4
c = 0.53
C_Dp = 0.020
Omega = 30.264 #rad/s
a = 343
CdS = 1.65
vi_hover = 15.06
V = np.arange(1,82,1)
#-------- Base calculations ----------#
sigma = N_blades*c/(pi*R)
V_tip = Omega*R/60
Mach_tip = V_tip/a
mu = V/(Omega*R)

#-------- Main calculations ----------#
def Powercurves(V):
    V_bar = V/vi_hover
    vi_bar = 1/V_bar
    Pi = k*W*vi_bar*sqrt(W/(2*rho*pi*R**2))
    Pp = (sigma*C_Dp)/8*rho*(Omega*R)**3*pi*R**2*(1+4.65*mu**2)
    Ppar = CdS*0.5*rho*V**3
    return Pi, Pp, Ppar

Pil = []
Ppl = []
Pparl = []

for i in range(len(V)):
    Pi, Pp, Ppar = Powercurves(V[i])
    Pil.append(Pi)
    Ppl.append(Pp)
    Pparl.append(Ppar)

# print(Pil)


#-------- Plotting ----------#
plt.plot(V,Pil)