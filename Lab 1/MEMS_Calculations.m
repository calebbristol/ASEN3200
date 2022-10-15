%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ASEN 3200 Lab 1 Calculations Code  %
%                                      %
%   Author: Caleb Bristol              %
%   Date: 09/15/21                     %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminary MATLAB

clc
clear
close all;

%% Read in Data
MEMS_MAN = load("2021_09_09_TH_UNIT05_MEMS_MAN");
MEMS_F2C1 = load("2021_09_09_TH_UNIT05_MEMS_F2C1");
MEMS_F2C5 = load("2021_09_09_TH_UNIT05_MEMS_F2C5");
MEMS_F5C5 = load("2021_09_09_TH_UNIT05_MEMS_F5C5");

MEMS_MAN(1,:) = [];
MEMS_F2C1(1,:) = [];
MEMS_F2C5(1,:) = [];
MEMS_F5C5(1,:) = [];

MEMS_MAN(:,1) = MEMS_MAN(:,1) - MEMS_MAN(1,1);
MEMS_F2C1(:,1) = MEMS_F2C1(:,1) - MEMS_F2C1(1,1);
MEMS_F2C5(:,1) = MEMS_F2C5(:,1) - MEMS_F2C5(1,1);
MEMS_F5C5(:,1) = MEMS_F5C5(:,1) - MEMS_F5C5(1,1);

%% Convert RPM for input rate into rad/s

MEMS_MAN(:,3) = (1/60) * (2*pi) * MEMS_MAN(:,3);
MEMS_F2C1(:,3) = (1/60) * (2*pi) * MEMS_F2C1(:,3);
MEMS_F2C5(:,3) = (1/60) * (2*pi) * MEMS_F2C5(:,3);
MEMS_F5C5(:,3) = (1/60) * (2*pi) * MEMS_F5C5(:,3);

%% Create Line of Best Fit for Gyro vs. Encoder

bf_F2C1 = polyfit(MEMS_F2C1(:,3),MEMS_F2C1(:,2),1);
bf_F2C5 = polyfit(MEMS_F2C5(:,3),MEMS_F2C5(:,2),1);
bf_F5C5 = polyfit(MEMS_F5C5(:,3),MEMS_F5C5(:,2),1);
bf_MAN = polyfit(MEMS_MAN(:,3),MEMS_MAN(:,2),1);

lbf_F2C1 = bf_F2C1(1) * MEMS_F2C1(:,3) + bf_F2C1(2);
lbf_F2C5 = bf_F2C5(1) * MEMS_F2C5(:,3) + bf_F2C5(2);
lbf_F5C5 = bf_F5C5(1) * MEMS_F5C5(:,3) + bf_F5C5(2);
lbf_MAN = bf_MAN(1) * MEMS_MAN(:,3) + bf_MAN(2);

%% Plot Gyro vs. Encoder
figure(1)

plot(MEMS_F2C1(:,3),MEMS_F2C1(:,2),'LineWidth',1.5); hold on
plot(MEMS_F2C1(:,3),lbf_F2C1,'LineWidth',1.5)
title("Gyro Output vs. Encoder Velocity (0.2 Hz, 0.1 A)")
xlabel("Encoder Velocity [rad/s]")
ylabel("Gyro Output [rad/s]")
legend("Data","Line of Best Fit")
set(gca,'fontsize',12)
hold off

figure(2)

plot(MEMS_F2C5(:,3),MEMS_F2C5(:,2),'LineWidth',1.5); hold on
plot(MEMS_F2C5(:,3),lbf_F2C5,'LineWidth',1.5)
title("Gyro Output vs. Encoder Velocity (0.2 Hz, 0.5 A)")
xlabel("Encoder Velocity [rad/s]")
ylabel("Gyro Output [rad/s]")
legend("Data","Line of Best Fit")
set(gca,'fontsize',12)
hold off

figure(3)

plot(MEMS_F5C5(:,3),MEMS_F5C5(:,2),'LineWidth',1.5); hold on
plot(MEMS_F5C5(:,3),lbf_F5C5,'LineWidth',1.5)
title("Gyro Output vs. Encoder Velocity (0.5 Hz, 0.5 A)")
xlabel("Encoder Velocity [rad/s]")
ylabel("Gyro Output [rad/s]")
legend("Data","Line of Best Fit")
set(gca,'fontsize',12)
hold off

figure(4)

plot(MEMS_MAN(:,3),MEMS_MAN(:,2),'LineWidth',1.5); hold on
plot(MEMS_MAN(:,3),lbf_MAN,'LineWidth',1.5)
title("Gyro Output vs. Encoder Velocity (Manual Control)")
xlabel("Encoder Velocity [rad/s]")
ylabel("Gyro Output [rad/s]")
legend("Data","Line of Best Fit")
set(gca,'fontsize',12)
hold off

%% Summate Delta T to create Theta 

delta_t = MEMS_MAN(2,1) - MEMS_MAN(1,1);

theta_F2C1 = cumsum(MEMS_F2C1(:,2) - bf_F2C1(2)) .* delta_t ./ bf_F2C1(1);
theta_F2C1r = cumsum(MEMS_F2C1(:,3)) .* delta_t;
theta_F2C5 = cumsum(MEMS_F2C5(:,2) - bf_F2C5(2)) .* delta_t ./ bf_F2C5(1);
theta_F2C5r = cumsum(MEMS_F2C5(:,3)) .* delta_t;
theta_F5C5 = cumsum(MEMS_F5C5(:,2) - bf_F5C5(2)) .* delta_t ./ bf_F5C5(1);
theta_F5C5r = cumsum(MEMS_F5C5(:,3)) .* delta_t;
theta_MAN = cumsum(MEMS_MAN(:,2) - bf_MAN(2)) .* delta_t ./ bf_MAN(1);
theta_MANr = cumsum(MEMS_MAN(:,3)) .* delta_t;

%% Plot Angular Position
figure(5)

plot(MEMS_F2C1(:,1),theta_F2C1,'LineWidth',1.5); hold on 
plot(MEMS_F2C1(:,1),theta_F2C1r,'--','LineWidth',1.5)
title("Angular Position vs. Time (0.2 Hz, 0.1 A)")
xlabel("Time [s]")
ylabel("Angular Position \theta [rad]")
legend("Fitted Gyro Data","Truth Data",'location','NW')
set(gca,'fontsize',12)
hold off

figure(6)

plot(MEMS_F2C5(:,1),theta_F2C5,'LineWidth',1.5); hold on 
plot(MEMS_F2C5(:,1),theta_F2C5r,'--','LineWidth',1.5)
title("Angular Position vs. Time (0.2 Hz, 0.5 A)")
xlabel("Time [s]")
ylabel("Angular Position \theta [rad]")
legend("Fitted Gyro Data","Truth Data",'location','NW')
set(gca,'fontsize',12)
hold off

figure(7)

plot(MEMS_F5C5(:,1),theta_F5C5,'LineWidth',1.5); hold on 
plot(MEMS_F5C5(:,1),theta_F5C5r,'--','LineWidth',1.5)
title("Angular Position vs. Time (0.5 Hz, 0.5 A)")
xlabel("Time [s]")
ylabel("Angular Position \theta [rad]")
legend("Fitted Gyro Data","Truth Data")
set(gca,'fontsize',12)
hold off

figure(8)

plot(MEMS_MAN(:,1),theta_MAN,'LineWidth',1.5); hold on 
plot(MEMS_MAN(:,1),theta_MANr,'--','LineWidth',1.5)
title("Angular Position vs. Time (Manual Control)")
xlabel("Time [s]")
ylabel("Angular Position \theta [rad]")
legend("Fitted Gyro Data","Truth Data",'location','NW')
set(gca,'fontsize',12)
hold off

%% Std and Mean

Std_K = std([bf_F2C1(1) bf_F2C5(1) bf_F5C5(1)]);
Mean_K = mean([bf_F2C1(1) bf_F2C5(1) bf_F5C5(1)]);

Std_B = std([bf_F2C1(2) bf_F2C5(2) bf_F5C5(2)]);
Mean_B = mean([bf_F2C1(2) bf_F2C5(2) bf_F5C5(2)]);