%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% AE4314 Assignment %%%%%%%%%%
%%% Yara Hinssen and Sybren Bootsma %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%

%Trim function 
function [a_1, theta0,lambda_i] = trim_con(V) 
    W = 93440.25 ; %N
    rho = 1.225 ;
    R = 7.315 ; %m
    Omega = 30.264 ; %rad/s
    CdS = 1.65 ;
    N_blades = 4 ;
    c = 0.53 ;
    Cl_alpha = 0.11*180/pi  ;%Fairfoiltools
    D = 0.5 * rho * CdS * V^2 ;
    T = sqrt(W^2 + D^2) ;
    sigma = N_blades*c/(pi*R) ;
    V_tip = Omega*R ;
    C_T = T / (rho * pi * R^2 * V_tip^2) ; 
    vi_0 = sqrt(W / (2*rho*pi*R^2)) ;
    lambda_i = vi_0 / V_tip ; %starting value for lambda_i
    mu_trim = V / V_tip ; 
    error = 1 ;
    epsilon = 0.0001;
    while error > epsilon
        C_T_glau = 2*lambda_i * sqrt((mu_trim * cos(D/W))^2 + (mu_trim * sin(D/W) + lambda_i)^2) ;
        lambda_i = lambda_i - 0.000001 ; 
        error = abs(C_T - C_T_glau) ;    
    end
    A = [(1 + 1.5*mu_trim^2) (-8/3 * mu_trim); -mu_trim (2/3 + mu_trim^2)] ;
    a_1 = ((-2*mu_trim^2 * D/W - 2*mu_trim*lambda_i) / det(A)) * 180/pi ;
    theta0 = ((4/sigma * C_T / Cl_alpha + mu_trim*D/W + lambda_i) / det(A)) * 180/pi ;
end
