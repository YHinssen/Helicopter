%SIMULATION OF A CYCLIC PITCH INPUT THETA_C=1DEG GIVEN FROM HOVER at
%0.5 sec<t<1 sec as in the lecture notes  
clear
%INITIAL DATA HELICOTPER
g = 9.81;	
Cl_alpha = 5.7; %NACA 0012
sigma = .075;	%blade solidity	
gamma = 6;
CdS = 1.5;
m = 2200;
rho = 1.225;
vtip = 200;
R = 7.32;
Myy = 10615;
h_rotor = 1;
Omega = vtip /(R);
A = pi*R^2;
tau = .1;		%time constant in dynamiCs inflow!!!
collect(1) = 6*pi/180;
longit(1) = 0*pi/180;

%INITIAL VALUES SIMULATION;
t0=0;
u0=0;
w0=0;
q0=0;
pitch0=0*pi/180;
x0=0;
labi0=sqrt(m*g/(A*2*rho))/vtip;%nondimentional inflow

t(1)=t0;
u(1)=u0;
w(1)=w0;
q(1)=q0;
pitch(1)=pitch0;
x(1)=x0;
labi(1)=labi0;
z(1)=0;

%START INTEGRATION FOR 40 SECONDS
aantal=400;
teind=40;
stap=(teind-t0)/aantal;

for i=1:aantal 

   if t(i)>=0.5 & t(i)<=1 longit(i)=1*pi/180;%longit cyclic input
   else longit(i)=0*pi/180;
   end
   
   %no LAW FOR COLLECTIVE
c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i));
h(i)=-z(i);
collect(i)=collect(1);    

%Defining the differential equations

%defining the nondimensional notations
qdiml(i)=q(i)/Omega;
vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip;
if u(i)==0 	if w(i)>0 	phi(i)=pi/2;
        else phi(i)=-pi/2;end
else
phi(i)=atan(w(i)/u(i));
end
if u(i)<0
phi(i)=phi(i)+pi;
end
alfc(i)=longit(i)-phi(i);

mu(i)=vdiml(i)*cos(alfc(i));
labc(i)=vdiml(i)*sin(alfc(i));

%a1 Flapping calculus
teller(i)=-16/gamma*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
a1(i)=teller(i)/(1-.5*mu(i)^2);

%the thrust coefficient
ctelem(i)=Cl_alpha*sigma/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));
%Thrust coefficient from Glauert
alfd(i)=alfc(i)-a1(i);
ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*...
sin(alfd(i))+labi(i))^2);

%Equations of motion
labidot(i)=ctelem(i); 
thrust(i)=labidot(i)*rho*vtip^2*A;
helling(i)=longit(i)-a1(i);
vv(i)=vdiml(i)*vtip; 		%it is 1/sqrt(u^2+w^2)

udot(i)=-g*sin(pitch(i))-CdS/m*.5*rho*u(i)*vv(i)+...
thrust(i)/m*sin(helling(i))-q(i)*w(i);

wdot(i)=g*cos(pitch(i))-CdS/m*.5*rho*w(i)*vv(i)-...
thrust(i)/m*cos(helling(i))+q(i)*u(i);

qdot(i)=-thrust(i)*h_rotor/Myy*sin(helling(i));

pitchdot(i)=q(i);

xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i));

zdot(i)=-c(i);

labidot(i)=(ctelem(i)-ctglau(i))/tau;

%Euler integration
u(i+1)=u(i)+stap*udot(i);
w(i+1)=w(i)+stap*wdot(i);
q(i+1)=q(i)+stap*qdot(i);
pitch(i+1)=pitch(i)+stap*pitchdot(i);
x(i+1)=x(i)+stap*xdot(i);
labi(i+1)=labi(i)+stap*labidot(i);
z(i+1)=z(i)+stap*zdot(i);
t(i+1)=t(i)+stap;
end;


figure(1)
plot(t(1:400),longit*180/pi),xlabel('t (s)'),ylabel('longit grd')
figure(2)
plot(t,u),xlabel('t (s)'),ylabel('u(m)')
figure(3)
plot(t,pitch*180/pi), xlabel('t (s)'), ylabel('pitch(deg)')
figure(4)
plot(t,x),xlabel('t (s)'),ylabel('x(m)')
figure(5)
plot(t,w),xlabel('t (s)'),ylabel('w(m)')
figure(6)
plot(t,q*180/pi),xlabel('t (s)'),ylabel('q(m)')
figure(7) 
plot(t,labi),xlabel('t (s)'),ylabel('labi(m)')
figure(8)
plot(t,-z),xlabel('t (s)'),ylabel('h(m)')

S1 = stepinfo(pitch*180/pi,t);
S2 = stepinfo(q*180/pi,t);




