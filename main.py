import math as m
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
V = np.arange(1,100,1)
#-------- Base calculations ----------#
sigma = N_blades*c/(m.pi*R)
V_tip = Omega*R/60
Mach_tip = V_tip/a

#-------- Main calculations ----------#
def Powercurves(V, vi_bar):
    mu = V/(Omega*R)
    V_bar = V/vi_hover
    Pi = k*W*vi_bar*m.sqrt(W/(2*rho*m.pi*R**2))
    Pp = (sigma*C_Dp)/8*rho*(Omega*R)**3*m.pi*R**2*(1+4.65*mu**2)
    Ppar = CdS*0.5*rho*V**3
    return Pi, Pp, Ppar

Pil = []
Ppl = []
Pparl = []
Ptotal = []

for i in range(len(V)):
    if V[i] <= 12:
        V_bar = V[i]/vi_hover
        vi_bar = m.sqrt(((-V_bar)**2 / 2) + m.sqrt(((V_bar)**4 / 4) + 1))
    else: 
        V_bar = V[i]/vi_hover
        vi_bar = 1/V_bar
    Pi, Pp, Ppar = Powercurves(V[i], vi_bar)
    Pil.append(Pi)
    Ppl.append(Pp)
    Pparl.append(Ppar)
    Ptotal.append(Pi+Pp+Ppar)

# print(Pil)


#-------- Plotting ----------#
plt.figure(1)
plt.plot(V,Pil, label="$P_i$")
plt.plot(V,Ppl, label="$P_p$")
plt.plot(V,Pparl, label="$P_{par}$")
plt.plot(V, Ptotal, label="$P_{tot}$")
plt.legend(loc="upper center")
plt.show