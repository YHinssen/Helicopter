%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AE4314 Assignment %%%%%%%%%%
%%% Yara Hinssen and Sybren Bootsma %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%
k = 1.15 ;
W = 93440.25 ; %N
rho = 1.225 ;
R = 7.315 ; %m
N_blades = 4 ;
c = 0.53 ;
C_Dp = 0.015 ;
Omega = 30.264 ; %rad/s
a = 343 ;
CdS = 1.65 ;
vi_hover = 15.06 ;
V = 1:0.01:100 ;
t = linspace(0,(2*pi/Omega),200) ;
g = 9.80665 ;

% ------ Constants Tail Rotor---------%
R_tr = 2.79/2  ;
c_tr = 0.833 * 0.3048 ;   %m 
N_tr = 4 ;
l_tr = 8.998 ;            %m
Omega_tr = 1530 / 2 / pi ;  %rad/s
sigma_tr = N_tr*c_tr/(pi*R_tr) ;

% ------ Constants Flapping-----------%
Cl_alpha = 0.11*180/pi  ;%Fairfoiltools
beta_0 = 0.2 ; %rad
theta_0 = 12.53/180*pi ; %rad
q = 20/180*pi ; %rotation rate
p = 10/180*pi ;

%-------- Base calculations ----------%
sigma = N_blades*c/(pi*R) ;
V_tip = Omega*R ;
Mach_tip = V_tip/a ;

%% Moment of Inertia calculations

m_rotor = 4 * 110 ; %estimated from Sikorsky UH-60 helicopter
Myy_rotor = 1/12 * m_rotor * (c^2 + (2*R)^2); 
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
frac = (10.993 - ((x_cg * W / 9.80665 + 4.967*m_rotor + 14.060*m_tailrotor)/m_body) )/6.757 ;

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

%% Flapping 

%%%%% Flapping

m_blade = 110 ; %based on Blackhawk blade weight
i_blade = 1/3*m_blade*R^2 ; 
inflow_r = vi_hover/(Omega*R);
lock_nr = (rho*Cl_alpha*c*R^4)/i_blade ;
Vflap = 20 ;
alpha_c = 12.83/180*pi ;
mu = Vflap*sin(alpha_c)/(Omega*R) ;
lam_c = Vflap*sin(alpha_c)/(Omega*R) ;

%Angular frame
a0 = lock_nr/8*(theta_0-(4/3)*inflow_r) ;
a1 = p/Omega -16/lock_nr*(q/Omega) ;
b1 = -q/Omega - 16*p/(lock_nr*Omega);

a0forward = lock_nr/8*(theta_0*(1+mu^2)-4/3*(lam_c+inflow_r)) ;
a1forward = (8/3*mu*theta_0-2*mu*(lam_c+inflow_r)-16/lock_nr*q/Omega)/(1-1/2*mu^2) ;
b1forward = (4/3*mu*a0forward-q/Omega)/(1+1/2*mu^2) ;

beta_par = lock_nr/8*(theta_0-(4*inflow_r)/3)*180/pi ;

beta = zeros(length(t),1) ;
beta_rot = zeros(length(t),1) ;
beta_homl = zeros(length(t),1) ;
beta_for = zeros(length(t),1) ;
for i = 1: length(t)
    beta_hom = (beta_0 * exp(-lock_nr/16*Omega*t(i))*(cos(Omega*sqrt(1-(lock_nr/16)^2)*t(i)) + (lock_nr/16)/sqrt(1-(lock_nr/16)^2) * sin(Omega*sqrt(1-(lock_nr/16)^2)*t(i)))*180/pi) ;
    beta_homl(i) = beta_hom ;
    beta(i) = (beta_hom + beta_par);
    
    beta_parrot = (a0 - a1 * cos(Omega * t(i)) - b1 * sin(Omega * t(i)))*180/pi ;
    beta_rot(i) = (beta_hom + beta_parrot) ;
        
    beta_parfor = (a0forward - a1forward * cos(Omega * t(i)) - b1forward * sin(Omega * t(i)))*180/pi ;
    beta_for(i) = (beta_hom +beta_parfor) ;
end

figure(1)
plot(t, beta, 'DisplayName', "Hover")
hold on
yline(beta_par, 'DisplayName', "Hover Particular")
hold on
plot(t, beta_homl, 'DisplayName', "Hover Homogeneous")
hold on
plot(t,beta_rot, 'DisplayName', "Pitch and roll")
hold on
plot(t,beta_for, 'DisplayName', "Forward flight")
hold off
ylabel("Beta [degrees]")
xlim([0,t(end)])
legend


%% Trim calculations

a1 = zeros(length(V),1) ;
theta_0 = zeros(length(V),1) ;
epsilon = 0.001;

vi_0 = sqrt(W / (2*rho*pi*R^2)) ;
lambda_i = vi_0 / V_tip ; %starting value for lambda_i
for i = 1: length(V)
    D = 0.5 * rho * CdS * V(i)^2 ;
    T = sqrt(W^2 + D^2) ;
    C_T = T / (rho * pi * R^2 * V_tip^2) ; 
    mu_trim = V(i) / V_tip ; 
    lambda_i = vi_0 / V_tip ; %starting value for lambda_i
    error = 1 ;
    while error > epsilon
        C_T_glau = 2*lambda_i * sqrt((mu_trim * cos(D/W))^2 + (mu_trim * sin(D/W) + lambda_i)^2) ;
        lambda_i = lambda_i - 0.000001 ; 
        error = abs(C_T - C_T_glau) ;    
    end
    
    A = [(1 + 1.5*mu_trim^2) (-8/3 * mu_trim); -mu_trim (2/3 + mu_trim^2)] ;
    a1(i) = ((-2*mu_trim^2 * D/W - 2*mu_trim*lambda_i) / det(A)) * 180/pi ;
    theta_0(i) = ((4/sigma * C_T / Cl_alpha + mu_trim*D/W + lambda_i) / det(A)) * 180/pi ;
end

% Plotting
figure(2)
plot(V, -a1, 'DisplayName', '\theta_c')
hold on
plot(V, theta_0, 'DisplayName', '\theta_0')
hold off
legend

%% Manoeuver Simulation

V1 = 90  * 0.514444 ;
V2 = 70  * 0.514444 ;
V3 = 90  * 0.514444 ;
V4 = 110 * 0.514444 ;
m = W/g ;
h = 2.045 ; %height of center of the diskplane above cg

%Store info
t = [] ;
V = [] ;
u = [] ;
w = [] ;
q = [] ;
theta_f = [] ;
x = [];
z = [];


%Set initial values
t0 = 0 ;
V(1) = V1 ;
t(1) = t0 ;
[a_1, theta0] = trim_con(V1) ;
x(1) = 0;

i = 1 ;
test = 1 ;
ref = 1 ;

%Set time step
dt = 0.1 ;
while test == 1
    i = i + 1 ;
    t(i) = t(i-1) + dt ;
    if ref == 1
        V_ref = V2;
    end
%     if ref = 2
%         V_ref = V3;
%     end 
%     if ref = 3
%         V_ref = V4;
%     end
       
    V(i) = sqrt(u(i-1)^2 + w(i-1)^2) ; 
    D = 0.5 * rho * CdS * V(i)^2 ;
    
    theta_f_ref = 5 * 180 / pi ;
    %-0.005 * (15 - x) + 0.02 * u(i-1) ;
    theta_c(i) = 0.2*(theta_f(i-1) - theta_f_ref) + 0.2*q(i-1);
    
    alpha_c = theta_c - arctan(w(i-1)/u(i-1)) ;
    mu = V(i) / (Omega*R) * cos(alpha_c) ;
    lambda_c = V / (Omega*R) * sin(alpha_c) ;
    
    lambda_i = vi_0 / V_tip ; %starting value for lambda_i
    error = 1 ;
    while error > epsilon
        a1 = (8/3 * mu * theta_0 - 2*mu * (lambda_c + lambda_i) - 16 / locknr * q(i-1)/Omega) / (1 - 0.5*mu^2) ; 
        C_T_BEM = 0.25 * Cl_alpha * sigma* (2/3*theta_0*(1+3/2*mu^2)-(lambda_c+lambda_i)) ;
        C_T_glau = 2*lambda_i * sqrt((mu_trim * cos(D/W))^2 + (mu_trim * sin(D/W) + lambda_i)^2) ;
        lambda_i = lambda_i - 0.000001 ; 
        error = abs(C_T_BEM - C_T_glau) ; 
    end
    
    T = C_T_glau * rho * (Omega * R)^2 * pi * R^2;
    
    udot = -g * sin(theta_f(i-1)) - D/m * u(i-1) / V(i-1) + T/m * sin(theta_c - a1) - q(i-1)*w(i-1);
    wdot = g * cos(theta_f(i-1)) - D/m * w(i-1) / V(i-1) + T/m * cos(theta_c - a1) + q(i-1)*u(i-1);
    qdot = -T/Myy_Total * h * sin(theta_c - a1);
    theta_fdot = q(i-1);
    xdot = u(i)*cos(theta_f(i))+w(i)*sin(theta_f(i));
    zdot = -u(i)*sin(theta_f(i))+w(i)*cos(theta_f(i));
    
    
    u(i)=u(i-1)+dt*udot;
    w(i)=w(i-1)+dt*wdot;
    q(i)=q(i-1)+dt*qdot;
    theta_f(i)=theta_f(i-1)+dt*theta_fdot;
    x(i)=x(i-1)+dt*xdot;
    z(i)=z(i-1)+dt*-zdot;
    
%     if reached
%         test = 0 ;
    test = 0 ;
end
