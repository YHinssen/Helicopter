%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AE4314 Assignment %%%%%%%%%%
%%% Yara Hinssen and Sybren Bootsma %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Phugoid calculation %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%
k = 1.15 ;
W = 93440.25 ; %N
m = W / 9.81 ;
rho = 1.225 ;
R = 7.315 ; %m
N_blades = 4 ;
c = 0.53 ;
C_Dp = 0.015 ;
Omega = 30.264 ; %rad/s
a = 343 ;
CdS = 1.65 ;
vi_hover = 15.06 ;
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
m_blade = 110 ; %based on Blackhawk blade weight
i_blade = 1/3*m_blade*R^2 ; 
lock_nr = (rho*Cl_alpha*c*R^4)/i_blade ;

%-------- Base calculations ----------%
sigma = N_blades*c/(pi*R) ;
V_tip = Omega*R ;
Mach_tip = V_tip/Cl_alpha ;

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

%% Phugoid calculations

V = 90 * 0.514444 ;
h_rotor = 2.045/R ; %height of center of the diskplane above cg
tau_c = 0 ; %level flight so climbing angle is zero
D = 0.5 * rho * CdS * V^2 ;
T = sqrt(D^2 + W^2) ;
C_T = T / (rho * pi * R^2 * V_tip^2);
t_c = T / (rho * pi * R^2 * V_tip^2 * sigma) ;
t_hat = m/(rho*sigma*pi*R^3 * Omega) ;
wc = W / (rho * pi * R^2 * V_tip^2 * sigma) ;
[a1, theta0, lambda_i] = trim_con(V) ;
% a1 = 2*mu*(4*theta0/3 + lambda) / (1 - 0.5*mu^2) ;
a1 = 0;
theta0 = 0 ;
theta_f0 = atan( - D/W);


mu = V/V_tip ;
delta = C_Dp  ; %blade profile drag coefficient
V_bar = V / vi_hover ;
V_hat = V / V_tip ;
nu_i_bar = sqrt((-(V_bar)^2 / 2) + sqrt(((V_bar)^4 / 4) + 1)) ; %

S_FP = 2.5; %CHECK EVEN MET PAINT %Fuselage equivalent plate area
d0 = S_FP / (sigma*pi*R^2) ;
hcd = 0.25*mu*delta ;
tcd = t_c ;
alpha_D = (-0.5 * mu^2 * d0 + hcd / tcd) ;
alpha_nf = alpha_D - a1;
lambda_D = mu * sin(alpha_D - lambda_i) ;
lambda = sin(alpha_nf - lambda_i) ;


l = (z_cg_rotor - z_cg)/R ;
b1 = a1 + hcd / t_c ;
a1_s = hcd / t_c ;
zq = 0 ; %done
a1 = 2*mu*(4*theta0/3 + lambda)/(1-0.5*mu^2);
theta0 = ((4 * tcd/Cl_alpha) - lambda_D * (1 - 0.5*mu^2) / (1+1.5*mu^2)) * 1.5*(1+1.5*mu^2)/(1-mu^2 + 2.25*mu^4) ;

dlambda_i_dmu = (2*mu*theta0 + alpha_nf - (4*t_c /(Cl_alpha * lambda_i))*V_bar*nu_i_bar^3) / (1 + (4/Cl_alpha)*(t_c/lambda_i)*(1 + nu_i_bar^4)) ;
dlambda_dmu = alpha_nf - dlambda_i_dmu ;
da1_dmu = a1 / mu - 2*mu/(1-0.5*mu^2) * dlambda_dmu ;
dtc_dmu = (2*mu*theta0 + alpha_D - a1 + V_bar * nu_i_bar^3 / (1 + nu_i_bar^4))/(4/Cl_alpha + (lambda_i/t_c)/(1 + nu_i_bar^4)) ;
dhcd_dmu = 0.25 * delta ;

da1_dw = 2*mu / ((1 - mu^2/2) * (1 + (Cl_alpha/4)*lambda_i/t_c + nu_i_bar^4)) ;
dtc_dw = Cl_alpha/4 /(1 + (Cl_alpha/4)*lambda_i/t_c + nu_i_bar^4);
dhcd_dw = Cl_alpha/4 * 1/(1 + (Cl_alpha/4)*lambda_i/t_c + nu_i_bar^4) * (0.5*a1 - mu*theta0 + mu*lambda_D/(1 - 0.5*mu^2)) ;
 
da1_dq = -16 / lock_nr /(1 - 0.5*mu^2) ;
dtc_dq = -zq ;
dhcd_dq = -4*Cl_alpha/(lock_nr * (1-0.5*mu^2)) * (0.5*lambda + mu*a1 - mu^2 * theta0) ;



xu = -t_c * da1_dmu - alpha_D * dtc_dmu - dhcd_dmu ;
xw = -t_c * da1_dw - alpha_D * dtc_dw - dhcd_dw ;
xq = -t_c * da1_dq - alpha_D * dtc_dq - dhcd_dq ;

zu = -dtc_dmu ;
zw = -dtc_dw  ;
zq = -dtc_dq  ;

mu_prime = - (l - h_rotor*a1_s) * dtc_dmu + h_rotor*(t_c * da1_dmu + dhcd_dmu) ;
mw_prime = - (l - h_rotor*a1_s) * dtc_dw + h_rotor*(t_c * da1_dw + dhcd_dw) ;
mq_prime = - (l - h_rotor*a1_s) * dtc_dq + h_rotor*(t_c * da1_dq + dhcd_dq) ;

mu_star = m/(rho*sigma*pi*R^3) ;
iB = Myy_total / (m * R^2) ;

m_u = mu_star / iB * mu_prime ;
m_w = mu_star / iB * mw_prime ;
m_q = -mq_prime / iB ;
% 
% dlambda_dw = alpha - dlambdai_dw ;
% dCT_dmu = Cl_alpha * sigma /4 * dlambda_dw ;
% dlambda0_dmu = alpha - 
% dlambdai0_dw = lambda_i0 / C_T * dCT_dmu 0 8 * lambda_i0^3 / C_T^2 * (mu + lambda0 * dlambda0_dmu)
% dlambdai_dw = (1 + k) * dlambdai0_dw ;
% m_wdot_prime = -V_inf * l_inf * a_inf * dlambdai_dw ;
m_wdot = 0 ;

% B matrix - control matrix
da1_dB1 = -mu * da1_dw   ;
dtc_dB1 = -mu * dtc_dw   ;
dhcd_dB1 = -mu * dhcd_dw ;

dlambdai_dtheta0 = (2/3*Cl_alpha*sigma*(1 + 1.5*mu^2))/(8 * mu + Cl_alpha * sigma) ;
da1_dtheta0 = 2*mu/(1 - 0.5*mu^2) * (4/3 - dlambdai_dtheta0) ;
dlambdaD_dtheta0 = mu*da1_dtheta0 - dlambdai_dtheta0 ;
dtc_dtheta0 = 4/3*Cl_alpha*mu*(1+1.5*mu^2)/(8*mu + Cl_alpha*sigma) ;
dhcd_dtheta0 = Cl_alpha/8 * ((a1 * dlambdaD_dtheta0 + lambda_D*da1_dtheta0) - 2*mu*(lambda_D + theta0 * dlambdaD_dtheta0)) ;

xB1 = dtc_dB1 * alpha_D + t_c * (1+ mu * da1_dw) - dhcd_dB1 ;
zB1 = -dtc_dB1 ;
mB1_prime = -(l-h_rotor*a1_s)*dtc_dB1 - t_c*h_rotor * (1 + mu * da1_dw) + h_rotor * dhcd_dB1 ;
mB1 = mu_star*mB1_prime/iB ; 

xtheta0 = -t_c* da1_dtheta0 - alpha_D*dtc_dtheta0 - dhcd_dtheta0 ; 
ztheta0 = -dtc_dtheta0 ;
mtheta0_prime = -(l - h_rotor * a1_s)* dtc_dtheta0 + t_c*h_rotor * da1_dtheta0 + h_rotor* dhcd_dtheta0 ;
mtheta0 = mu_star*mtheta0_prime/iB ;



Ktheta = 0.75 ; %From Code.m --> PID controller
Kq = 0.75 ;

A = [xu xw (-C_T*cos(tau_c)) xq; zu zw (-C_T*sin(tau_c)) (mu+zq); 0 0 0 1; ...
    (m_u+zu*m_wdot) (m_w + zw*m_wdot) (-m_wdot*C_T*sin(tau_c)) (m_q + m_wdot*(V_hat + zq))] ;
B = [xB1 zB1 0 (mB1 + m_wdot*zB1);
    xtheta0 ztheta0 0  (mtheta0 + m_wdot*ztheta0)] ;
K = [0 0 Ktheta Kq/t_hat; 0 0 0 0] ;

eigen_nocontrol = eig(A) ;
eigen_control = eig(A-B.' * K) ;

% xq_prime = 
% zq_prime = 
% xthetaf = 
% zthetaf = 
% mthetaf_prime = 

Xu = xu * rho * sigma * pi * R^3 * Omega ;
% Xw = xw * rho * sigma * pi * R^3 * Omega ;
% Xq = xq_prime * rho * sigma * pi * R^4 * Omega ;
% Xthetaf = xthetaf * rho * sigma * pi * R^4 * Omega^2 ;
% XB1 = xB1 * rho * sigma * pi * R^4 * Omega^2 ;
% Xtheta0 = xtheta0 * rho * sigma * pi * R^4 * Omega^2 ;

% Zu = zu * rho * sigma * pi * R^3 * Omega ;
% Zw = zw * rho * sigma * pi * R^3 * Omega ;
% Zq = zq_prime * rho * sigma * pi * R^4 * Omega ;
% Zthetaf = zthetaf * rho * sigma * pi * R^4 * Omega^2 ;
% ZB1 = zB1 * rho * sigma * pi * R^4 * Omega^2 ;
% Ztheta0 = ztheta0 * rho * sigma * pi * R^4 * Omega^2 ;

Mu = mu_prime * rho * sigma * pi * R^4 * Omega ;
% Mw = mw_prime * rho * sigma * pi * R^4 * Omega ;
Mq = mq_prime * rho * sigma * pi * R^5 * Omega ;
% Mthetaf = mthetaf_prime * rho * sigma * pi * R^5 * Omega^2 ;
% MB1 = mtheta0_prime * rho * sigma * pi * R^5 * Omega^2 ;
% Mtheta0 = mtheta0_prime * rho * sigma * pi * R^5 * Omega^2 ;

syms eigen
eqn = eigen^3 - (Xu/m + Mq/Myy_total) * eigen^2 + Xu/m*Mq/Myy_total * eigen + g*Mu/Myy_total == 0;
S = solve(eqn, eigen, 'MaxDegree', 3) ;

eig1 = simplify(S(1))
eig2 = simplify(S(2))
eig3 = simplify(S(3))

eig1_actual = (582538391040597808632503915665530953184334196401543 - 211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2))^(1/3)/216172782113783808 + (211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2) + 582538391040597808632503915665530953184334196401543)^(1/3)/216172782113783808 - 6893244655668713/216172782113783808 ;
eig2_actual = (3^(1/2)*(211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2) + 582538391040597808632503915665530953184334196401543)^(1/3)*1i)/432345564227567616 - (211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2) + 582538391040597808632503915665530953184334196401543)^(1/3)/432345564227567616 - ((582538391040597808632503915665530953184334196401543 - 211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2))^(1/3)*(1 + 3^(1/2)*1i))/432345564227567616 - 6893244655668713/216172782113783808 ;
eig3_actual = - (3^(1/2)*(211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2) + 582538391040597808632503915665530953184334196401543)^(1/3)*1i)/432345564227567616 - (211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2) + 582538391040597808632503915665530953184334196401543)^(1/3)/432345564227567616 + ((582538391040597808632503915665530953184334196401543 - 211106232532992*7614595609951795871339376229364821581342807978608078017583117661021705669^(1/2))^(1/3)*(- 1 + 3^(1/2)*1i))/432345564227567616 - 6893244655668713/216172782113783808 ;

omega_0 = imag(eig2_actual) ;
nu = - real(eig2_actual) ;
omega_n = sqrt(nu^2 + omega_0^2) ;
zeta = nu/omega_n ;
