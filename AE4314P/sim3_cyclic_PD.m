%SIMULATION OF A CYCLIC PITCH INPUT THETA_C=1 DEG GIVEN FROM HOVER 
%0.5 SEC<T<1 SEC. Now from the 15th second a PD controller becomes active 
clear
%INITIAL DATA HELICOPTER
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
tau = .1;		%time constant in dynamics inflow!!!
collect(1) = 6*pi/180;
longit(1) = 0*pi/180;



DL = m*g/A ;

DL_tab = [DL*0.7 DL*0.8 DL*0.9 DL DL*1.1 DL*1.2 DL*1.3; ...
    DL DL DL DL DL DL DL; DL DL DL DL DL DL DL; ...
    DL DL DL DL DL DL DL] ;

vtip_tab = [vtip vtip vtip vtip vtip vtip vtip; ...
    vtip*0.7 vtip*0.8 vtip*0.9 vtip vtip*1.1 vtip*1.2 vtip*1.3; ...
    vtip vtip vtip vtip vtip vtip vtip; vtip vtip vtip vtip vtip vtip vtip];

sigma_tab = [sigma sigma sigma sigma sigma sigma sigma; ...
    sigma sigma sigma sigma sigma sigma sigma; ...
    sigma*0.7 sigma*0.8 sigma*0.9 sigma sigma*1.1 sigma*1.2 sigma*1.3; ...
    sigma sigma sigma sigma sigma sigma sigma] ;

gamma_tab = [gamma gamma gamma gamma gamma gamma gamma; 
    gamma gamma gamma gamma gamma gamma gamma; ...
    gamma gamma gamma gamma gamma gamma gamma;
    gamma*0.7 gamma*0.8 gamma*0.9 gamma gamma*1.1 gamma*1.2 gamma*1.3] ;

delta_theta_pk_tab = zeros(4,7);
delta_theta_min_tab = zeros(4,7);
q_pk_tab = zeros(4,7);
thetapk_qpk_tab = zeros(4,7);


rows = size(q_pk_tab);
columns = size(q_pk_tab);


for j=1: rows(1)
    
    for k=1: columns(2)

    %initial values;
    t0=0;
    u0=0;
    w0=0;
    q0=0;
    pitch0=0*pi/180;
    x0=0;
    labi0=sqrt(DL_tab(j,k)/(2*rho))/vtip_tab(j,k);

    t(1)=t0;
    u(1)=u0;
    w(1)=w0;
    q(1)=q0;
    pitch(1)=pitch0;
    x(1)=x0;
    labi(1)=labi0;
    z(1)=0;

    %INTEGRATION 
    aantal=800;
    teind=80;
    stap=(teind-t0)/aantal;

    for i=1:aantal 
       if t(i)>=0.5 & t(i)<=1 longit(i)=1*pi/180;
       else longit(i)=0*pi/180;
       end

       if t(i)>=15 longitgrd(i)=.2*pitch(i)*180/pi+.2*q(i)*180/pi;%PD in deg
           longit(i)=longitgrd(i)*pi/180;	%in rad
       end    
       %longit(i)=longitgrd(i)*pi/180;	%in rad

    %NO LAW FOR COLLECTIVE

    c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i));
    h(i)=-z(i);
    collect(i)=collect(1);

    %Defining the differential equations

    %defining the nondimensional notations
    qdiml(i)=q(i)/Omega;
    vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip_tab(j,k);
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

    %a1 Flapping calculi
    teller(i)=-16/gamma_tab(j,k)*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
    a1(i)=teller(i)/(1-.5*mu(i)^2);

    %the thrust coefficient
    ctelem(i)=Cl_alpha*sigma_tab(j,k)/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));
    %Thrust coefficient from Glauert
    alfd(i)=alfc(i)-a1(i);
    ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*...
    sin(alfd(i))+labi(i))^2);

    %Equations of motion
    labidot(i)=ctelem(i); 
    thrust(i)=labidot(i)*rho*vtip_tab(j,k)^2*A;
    helling(i)=longit(i)-a1(i);
    vv(i)=vdiml(i)*vtip_tab(j,k); 		%it is 1/sqrt(u^2+w^2)

    udot(i)=-g*sin(pitch(i))-CdS/m*.5*rho*u(i)*vv(i)+...
    thrust(i)/m*sin(helling(i))-q(i)*w(i);

    wdot(i)=g*cos(pitch(i))-CdS/m*.5*rho*w(i)*vv(i)-...
    thrust(i)/m*cos(helling(i))+q(i)*u(i);

    qdot(i)=-thrust(i)*h_rotor/Myy*sin(helling(i));

    pitchdot(i)=q(i);

    xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i));

    zdot(i)=-c(i);

    labidot(i)=(ctelem(i)-ctglau(i))/tau;
    %corrdot(i)=uwens-u(i);
    %corrcdot(i)=cwens(i)-c(i);

    u(i+1)=u(i)+stap*udot(i);

    w(i+1)=w(i)+stap*wdot(i);
    q(i+1)=q(i)+stap*qdot(i);
    pitch(i+1)=pitch(i)+stap*pitchdot(i);
    x(i+1)=x(i)+stap*xdot(i);
    labi(i+1)=labi(i)+stap*labidot(i);
    z(i+1)=z(i)+stap*zdot(i);
    t(i+1)=t(i)+stap;
    end;

%     figure(1)
%     plot(t,u),xlabel('t (s)'),ylabel('u(m)');
    figure(2)
    plot(t,pitch*180/pi),xlabel('t (s)'),ylabel('pitch(deg)')
%     figure(3)
%     plot(t,x),xlabel('t (s)'),ylabel('x(m)')
%     figure(4)
%     plot(t,w),xlabel('t (s)'),ylabel('w(m)')
    figure(5)
    plot(t,q*180/pi),xlabel('t (s)'),ylabel('q(m)') 
%     figure(6)
%     plot(t,labi),xlabel('t (s)'),ylabel('labi(m)')
%     figure(7)
%     plot(t,-z),xlabel('t (s)'),ylabel('h(m)')
%     figure(8)
%     plot(t(1:800),longit*180/pi),xlabel('t (s)'),ylabel('longit grd')

    S1 = stepinfo(pitch*180/pi,t);
    S2 = stepinfo(q*180/pi,t);
    
    q_pk_tab(j,k) = S2.SettlingMax ;
    delta_theta_pk_tab(j,k) = S1.SettlingMax ;
    thetapk_qpk_tab(j,k) = S1.SettlingMax / S2.SettlingMax ;
    
    q_test = -100;
    for s=1: length(q)
        if q(s)*180/pi == q_pk_tab(j,k)
            q_test = q(s);
        end
        if q(s) < 0.9*q_test
            delta_theta_min_tab(j,k) = pitch(s)*180/pi;
            break
        end
    end
    
%     disp(q_pk_tab(j,k))
%     disp(q_test*180/pi)
%     disp(delta_theta_pk_tab(j,k))
%     disp(delta_theta_min_tab(j,k))
%     pause 
    end
end
