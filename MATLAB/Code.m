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
Omega_tr = 1530 / 2 / pi ;  %rad/s
sigma_tr = N_tr*c_tr/(pi*R_tr) ;

% Constants Flapping
Cl_alpha = 0.10867 ;
beta_0 = 0.2 ;
theta_0 = 12.53/180*pi ;

% Base calculations
sigma = N_blades*c/(pi*R) ;
V_tip = Omega*R ;
Mach_tip = V_tip/a ;

%% Moment of Inertia calculations

m_rotor = 4 * 110 ; %estimated from Sikorsky UH-60 helicopter
Myy_rotor = 0; %Negligible based on the small r from its own cg
x_cg_rotor = 4.967 ;
z_cg_rotor = 3.579 ; 

m_tailrotorblade = 0.58 * 14.594 ; %kg
m_tailrotor =  4 * 0.58 * 14.594 ; %kg
blade_cg = 1.6787 * 0.3048 * cos(27.5/180) ; %m
slug_ft2_to_kg_m2 = 1.35581795 ; %Conversion from slug ft^2 to kg m^2
Myy_tailrotor = 4 * (1.009 * slug_ft2_to_kg_m2 + blade_cg * m_tailrotorblade);
x_cg_tailrotor = 14.060 ;
z_cg_tailrotor = 2.885 ;

m_body = W / 9.80665 - m_rotor - m_tailrotor ;

%Determine fraction of body weight for body 1 and body 2
x_cg = 5.1562 ; %m
frac = (10.993 - ((cg_x * W / 9.80665 + 4.967*m_rotor + 14.060*m_tailrotor)/m_body) )/6.757 ;

m_body1 = m_body * frac ;
m_body2 = m_body * (1-frac) ;
Myy_body1 = 1/12 * (2.794^2 + 7.9249^2) * m_body1 ;
Myy_body2 = 1/12 * (1.716^2 + 6.5006^2) * m_body2 ;
x_cg_body1 = 4.236 ;
x_cg_body2 = 10.993 ;
z_cg_body1 = 1.607 ;
z_cg_body2 = 1.059 ;

z_cg = (m_rotor * z_cg_rotor + m_tailrotor * z_cg_tailrotor + m_body1 * z_cg_body1 + m_body2 * z_cg_body2) / (W / 9.80665);

r_rotor2 = (x_cg_rotor - x_cg)^2 + (z_cg_rotor - z_cg)^2 ;
r_tailrotor2 = (x_cg_tailrotor - x_cg)^2 + (z_cg_tailrotor - z_cg)^2 ;
r_body1_2 = (x_cg_body1 - x_cg)^2 + (z_cg_body1 - z_cg)^2 ;
r_body2_2 = (x_cg_body2 - x_cg)^2 + (z_cg_body2 - z_cg)^2 ;

Myy_total = Myy_rotor + m_rotor * r_rotor2 + Myy_body1 + m_body1 * r_body1_2 + Myy_body2 + m_body2 * r_body2_2 + Myy_tailrotor + m_tailrotor * r_tailrotor2 ;