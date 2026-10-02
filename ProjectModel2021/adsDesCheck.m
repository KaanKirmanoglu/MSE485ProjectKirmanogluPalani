clc
clear
close all

set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',30);
r_ads = 1.6;
r_des = 2;

t = 0:0.01:2.5;
theta_sim = readmatrix('surCov.txt');
theta_sim2 = readmatrix('surfCov3.txt');
theta = (r_ads/(r_ads+r_des))*(1-exp(-(r_ads+r_des)*t));
plot(t,theta,'-x','LineWidth',2)
hold on
%plot(t,theta_sim,'-o')
plot(t,theta_sim,'-o','LineWidth',2)
ylabel('\theta')
xlabel('time (s)')
legend('Analytical','Simulation')
title('r_{A} = 1.6, r_{D} = 2')
kb =1.38064852E-23;
N_av = 6.0221409E23; %Avogadro number
mO = 16/(1000*N_av);
mCO = 28/(1000*N_av);

