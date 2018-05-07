% Author : Sanjuksha Nirgude
%% Passivity based Adaptive Control 
clc;
clear all;
close all;
 
global torque
torque=[];

% Initial State
x0=[0.4,0.4,8,5,2.5];
tf=400; % Time Steps

options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4, 1e-4]);
%% IMPLEMENT THE CONTROLLER

[T,X] = ode45(@(t,x)Adaptive(t,x),[0 tf],x0,options);

%% Plot 
figure('Name','theta');
plot(T, X(:,1),'r-');
hold on
title('Theta under Adaptive Control');
plot(T,-sin(T),'b-');



figure('Name','Error in Inertia');
plot(T, X(:,3),'r-');
hold on
title('Error in Inertia in Adaptive Control');
plot(T,7.5*ones(size(T,1),1),'b-');

figure('Name','Error in Gravity ');
plot(T, X(:,4),'r-');
hold on
title('Error in Gravity in Adaptive Control');
plot(T,6*ones(size(T,1),1),'b-');

figure('Name','Error in Force');
plot(T, X(:,5),'r-');
hold on
title('Error in Force in Adaptive Control');
plot(T,1.5*ones(size(T,1),1),'b-');

figure('Name','Input_Adaptive control');
plot(T, torque(1,1:size(T,1)),'r-' );
title('Torque');

%% adaptive control

function [dx]=Adaptive(t,x) 
    % q q_dot  I_bar m_bar*g*d_bar fv_bar
    
    %x =[theta,dtheta, dI , dmgd , dfv];
    I = 7.5 ; mgd = 6; fv = 1.5;
    I_bar = x(3);
    mgd_bar = x(4);
    fv_bar = x(5);
    theta =x(1);
    dtheta =x(2);
    theta_d =[-sin(t)];
    dtheta_d =[-cos(t)];
    ddtheta_d =[sin(t)];
    M = I;
    C = fv;
    N = mgd;

    invM = inv(M);
    invMC = invM*C;

    M_bar = I_bar;
    C_bar = fv_bar;
    N_bar = mgd_bar;


    % Control Law
    e = theta -theta_d;
    de = dtheta -dtheta_d;
    
    Kv =475*eye(1);
    lamda = 1.9*eye(1);
    H =0.02*eye(3);
    
    a = ddtheta_d - lamda*de;
    v = dtheta_d - lamda*e;
    r = de + lamda*e;
    
    tau = M_bar * a + C_bar * v + N_bar - Kv * r;

    ddtheta = invM*tau - invM*N -invM*C*x(2);

    Y = [a, sin(theta), v];
    global torque
    torque=[torque,tau];
    
    dx=zeros(5,1);
    dx(1)= x(2); % dtheta
    dx(2)= ddtheta;
    dx(3:5)= -inv(H)*transpose(Y)*r;
end