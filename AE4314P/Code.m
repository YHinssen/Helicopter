%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AE4314P Helicopter Practical %%%%%
%%%% Yara Hinssen and Sybren Bootsma %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL DATA HELICOTPER
g = 9.81;	
Cl_alpha = 5.7; %NACA 0012
sigma = .075;	%blade solidity	
gamma = 6;
CdS = 1.5; % C_D * S
m = 2200; %kg
rho = 1.225; %Density kg/m3
vtip = 200;  %Tip speed m/s
R = 7.32;     %Blade radius
Myy = 10615;     %Mass moment of Inertia
h_rotor = 1;
Omega = vtip /(R);  %Rotational velocity
A = pi*R^2;
tau = .1;		%time constant in dynamics inflow
collect(1) = 6*pi/180;
longit(1) = 0*pi/180;
DL = m*g/A ;

