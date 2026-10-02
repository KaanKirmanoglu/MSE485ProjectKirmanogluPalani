clc
clear 
close all

temp = 1000:100:2000;
ssTotal= zeros(length(temp),1);
ssO  = zeros(length(temp),1);
ssCO =  zeros(length(temp),1);
oxiProb = zeros(length(temp),1);

for t=1:length(temp)
    filename1 = "surfCovTot" + string(temp(t)) + ".txt";
    surf_Covtotal = readmatrix(filename1);
    ssTotal(t) = mean(surf_Covtotal(3000:4000));
    filename2 = "surfCovO" + string(temp(t)) + ".txt";
    surf_CovO = readmatrix(filename2);
    ssO(t) = mean(surf_CovO(3000:4000));
    filename3 = "surfCovCO" + string(temp(t)) + ".txt";
    surf_CovCO = readmatrix(filename3);
    ssCO(t) = mean(surf_CovCO(3000:4000));
    filename4 = "carbonFlux" +string(temp(t)) + ".txt";
    carbonflux = readmatrix(filename4);
    oxiProb(t) = mean(carbonflux(3000:4000))/40; 
end

plot(temp,ssTotal,'x','LineWidth',2)
hold on
plot(temp,ssO,'o','LineWidth',2)
plot(temp,ssCO,'d','LineWidth',2)

figure

plot(temp,oxiProb,'x','LineWidth',2)




