clc
clear
close all

surfCovtotal = readmatrix('surfCovTot1200.txt');
surfO = readmatrix('surfCovO1200.txt');
surfCO = readmatrix('surfCovCO1200.txt');
t_arr = 0:1:length(surfCO)-1;

plot(t_arr,surfCovtotal,'-x')
hold on
plot(t_arr,surfO,'-d')
plot(t_arr,surfCO,'-o')

legend('tot','O','CO')