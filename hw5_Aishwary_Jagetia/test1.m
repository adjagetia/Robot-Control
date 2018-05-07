%% testing algorithm with a set of initial and final states.
clc
clear all
close all
%Example to demostrate the robust control with planar 2d arm.
theta10=-0.5;
dtheta10 =0; 
theta1f = 0.8;
dtheta1f=0;
tf=30;

% plan a trajectory to reach target postion given by theta1f, dot theta1f,
% theta2f, dot theta2f.
theta20=-1;
dtheta20= 0.1; 
theta2f = 0.5;
dtheta2f=0;

robustControl(theta10,theta20,dtheta10, dtheta20,theta1f, theta2f,dtheta1f,dtheta2f,tf)%plan_control(theta10,theta20,dtheta10, dtheta20,theta1f, theta2f,dtheta1f,dtheta2f,tf)
