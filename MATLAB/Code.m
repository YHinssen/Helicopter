%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AE4314 Assignment %%%%%%%%%%
%%% Yara Hinssen and Sybren Bootsma %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%
k = 1.15 ;
W = 93440.25 ;%N
rho = 1.225 ;
R = 7.315 ; %m
N_blades = 4 ;
c = 0.53 ;
C_Dp = 0.015 ;
Omega = 30.264 ; %rad/s
a = 343 ;
CdS = 1.65 ;
vi_hover = 15.06 ;
%V = np.arange(1,100,0.1)
%t = np.linspace(0,(6*m.pi/Omega),200)

% Constants Tail Rotor
R_tr = 2.79/2 ;
c_tr = 0.833 * 0.3048 ; %m 
N_tr = 4 ;
l_tr = 8.998 ;          %m
Omega_tr = 1530 / 2 / m.pi ;  %rad/s
sigma_tr = N_tr*c_tr/(m.pi*R_tr) ;

% Constants Flapping
Cl_alpha = 0.10867 ;
beta_0 = 0.2 ;
theta_0 = 12.53/180*m.pi ;

% Base calculations
sigma = N_blades*c/(m.pi*R) ;
V_tip = Omega*R ;
Mach_tip = V_tip/a ;

%% Moment of Inertia calculations

m_rotor = 4 * ((50*3.8 + 10*6.3)*10**(-2) * 4.45 / 0.27 / 9.81) ; %NASA doc -> [kg]
Myy_rotor = 0; %Negligible based on the small r from its own cg
z_cg_rotor = 

m_body = 
Myy_body = 
z_cg_body =

m_tailrotorblade = 0.58 * 14.594 ; %kg
m_tailrotor =  4 * 0.58 * 14.594 ; %kg
blade_cg = 1.6787 * 0.3048 * cos(27.5/180) ; %m
slug_ft2_to_kg_m2 = 1.35581795 ; %Conversion from slug ft^2 to kg m^2
Myy_tailrotor = 4 * (1.009 * slug_ft2_to_kg_m2 + blade_cg * m_tailrotorblade);
z_cg_tailrotor = 

Myy_total = Myy_rotor + m_rotor * z_cg_rotor^2 + Myy_body + m_body * z_cg_body^2 + Myy_tailrotor + m_tailrotor * z_cg_tailrotor^2 ;