clc
clear all
close all
A = [ 0 1 ; 0 0 ];
B = [ 0 ; 1 ];
K = [1 1.4];
temp_x=[];
temp_x_d=[];
x_0=[-3;1];
dt = 0.01;
for t=0:0.01:5
x1_d = - (0.0560*t^3) + (0.52*t^2) - 5;
x2_d = (1.04*t) - (0.1680*t^2);
x_d  = [x1_d; x2_d];    
u_d = 1.04 - (0.3360*t);
u_t=-K*(x_0-x_d)+u_d;
x_dot=A*x_0+B*u_t;
temp_x=[temp_x,x_0];
temp_x_d=[temp_x_d,x_d];
x_0=x_dot*dt+x_0;
end
figure

plot(temp_x(1,:),temp_x(2,:),'b-');
hold on
plot(temp_x_d(1,:), temp_x_d(2,:), 'r-');
legend('Desired','Actual')
title('Trajectory tracking controller');
xlabel('X-axis')
ylabel('Y-axis')
hold off

