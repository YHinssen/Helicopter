from math import *
import numpy as np
import matplotlib as plt

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
V = 0
#-------- Base calculations ----------#
sigma = N_blades*c/(pi*R)
V_tip = Omega*R/60
Mach_tip = V_tip/a
mu = V/(Omega*R)

#-------- Main calculations ----------#
def Powercurves(V):
    for i in range(0,84):
        V = V + 1
        V_bar = V/vi_hover
        vi_bar = 1/V_bar
        Pi = k*W*vi_bar*sqrt(W/(2*rho*pi*R**2))
        Pp = (sigma*C_Dp)/8*rho*(Omega*R)**3*pi*R**2*(1+4.65*mu**2)
        Ppar = CdS*0.5*rho*V**3
        Pil.append(Pi)
        Ppl.append(Pp)
        Pparl.append(Ppar)
    return Pil, Ppl, Pparl

Pil = []
Ppl = []
Pparl = []

Pil, Ppl, Pparl = Powercurves(V)

print(Pil)


#-------- Plotting ----------#
plt.plot(Pil,V)