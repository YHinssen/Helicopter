import math as m
import numpy as np
import matplotlib.pyplot as plt

# ----------- Constants --------------#
k = 1.15
W = 93440.25 #N
rho = 1.225
R = 7.315 #m
N_blades = 4
c = 0.53
C_Dp = 0.015
Omega = 30.264 #rad/s
a = 343
CdS = 1.65
vi_hover = 15.06
V = np.arange(1,100,0.1)

# ------ Constants Tail Rotor---------#
R_tr = 2.79/2
c_tr = 0.833 * 0.3048   #m 
N_tr = 4
l_tr = 8.998            #m
Omega_tr = 1530 / 2 / m.pi  #rad/s
sigma_tr = N_tr*c_tr/(m.pi*R_tr)


#-------- Base calculations ----------#
sigma = N_blades*c/(m.pi*R)
V_tip = Omega*R/60
Mach_tip = V_tip/a

#-------- Main calculations ----------#
def Powercurves(V):
    V_bar = V/vi_hover
    vi_bar = m.sqrt((-(V_bar)**2 / 2) + m.sqrt(((V_bar)**4 / 4) + 1))
    mu = V/(Omega*R)
    Pi = k*W*vi_bar*m.sqrt(W/(2*rho*m.pi*R**2))
    Pp = (sigma*C_Dp)/8 * rho * ((Omega*R)**3) * m.pi * (R**2) * (1+4.65*(mu**2))
    Ppar = CdS*0.5*rho*V**3
    
    
    T_tr = (Ppar + Pp + Pi)/(Omega * l_tr)
    mu_tr = V/(Omega_tr*R_tr)
    vi_bar_tr = m.sqrt((-(V_bar)**2 / 2) + m.sqrt(((V_bar)**4 / 4) + 1))
    Ptr = k*T_tr*vi_bar_tr*m.sqrt(T_tr/(2*rho*m.pi*R_tr**2)) + (sigma_tr*C_Dp)/8 * rho * ((Omega_tr*R_tr)**3) * m.pi * (R_tr**2) * (1+4.65*(mu_tr**2))
    return Pi, Pp, Ppar, Ptr

Pil = []
Ppl = []
Pparl = []
Ptrl = []
Ptotal = []

for i in range(len(V)):
    Pi, Pp, Ppar, Ptr = Powercurves(V[i])
    Pil.append(Pi/1000)
    Ppl.append(Pp/1000)
    Pparl.append(Ppar/1000)
    Ptrl.append(Ptr/1000)
    Ptotal.append((Pi+Pp+Ppar+Ptr)/1000)
    

#-------- Plotting ----------#
fig = plt.figure(1, figsize=(8,6))

ax = plt.subplot(111)
ax.plot(V,Pil, label="$P_i$")
ax.plot(V,Ppl, label="$P_p$")
ax.plot(V,Pparl, label="$P_{par}$")
ax.plot(V,Ptrl, label="$P_{tr}$")
ax.plot(V, Ptotal, label="$P_{tot}$")
plt.xlabel("V [m/s]")
plt.ylabel("P [kW]")
plt.xlim([0,100])
plt.ylim([0,2500])
plt.grid(True)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5)
plt.title("Power curve of the Boeing AH-64D Apache")
plt.show