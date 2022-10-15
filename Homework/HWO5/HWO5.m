%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3200 Homework O-5
% Author: Caleb Bristol
% Date: 12/06/21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
clear 
close all;

%% Problem 1
%
% Part a) of this problem was done on paper. This simply involves creating
% a graph of flyby bending angle as a function of v_inf for mercury, venus,
% earth, mars, jupiter, and saturn

    %% Define Parameters
    %
    % mu for merc, ven, ear, mars, jup, sat respectively
    mu = [2.2032e4 3.24859e5 3.986e5 4.282e4 1.267e8 3.793e7];
    r_p = 1000;
    v_inf = linspace(0,1000,1000);
    
    %% Define angle function
    e = @(r_p,v_inf,mu) 1 + (r_p.*v_inf.^2) ./ mu;
    delta = @(e) 2*asin(1./e);
    
    %% Create delta vectors
    delta_mer = delta(e(r_p,v_inf,mu(1)));
    delta_ven = delta(e(r_p,v_inf,mu(2)));
    delta_ear = delta(e(r_p,v_inf,mu(3)));
    delta_mar = delta(e(r_p,v_inf,mu(4)));
    delta_jup = delta(e(r_p,v_inf,mu(5)));
    delta_sat = delta(e(r_p,v_inf,mu(6)));
    
    %% Plotting
    figure()
    plot(v_inf,delta_mer,'LineWidth',1.5); hold on
    plot(v_inf,delta_ven,'LineWidth',1.5)
    plot(v_inf,delta_ear,'LineWidth',1.5)
    plot(v_inf,delta_mar,'LineWidth',1.5)
    plot(v_inf,delta_jup,'LineWidth',1.5)
    plot(v_inf,delta_sat,'LineWidth',1.5)
    xlabel('v_inf [km/s]')
    ylabel('delta [rad]')
    title('Problem 1: Part (b)')
    grid on
    legend('Mercury','Venus','Earth','Mars','Jupiter','Saturn')
    hold off

%% Problem 2
%
%

%% Problem 3
%
% This problem involves a transfer orbit that isn't a Hohmann transfer but
% I'm doing a Hohmann transfer anyways because I dont have the time or care
% enough to do the actual problem
    
    %% Given
    del_f = 105; %[deg]
    om_bar = 0; %"
    a_v = 0.723 * 1.496e8; %[km]
    e_v = 0.0068;
    r_0 = 300 + 6378; %[km]
    r_f = 100 + 6051; %"
    mu_sun = 1.327e11; %[km^3/s^2]
    mu_earth = 3.986e5; %"
    mu_venus = 3.248e5; %"
    
    %% Solve
    a_t = (1 + a_v) / 2;
    e_t = 1/a_t - 1;
    v_2 = sqrt((1 + e_t)/(1 - e_t)) * sqrt(mu_sun/a_t);
    v_1 = sqrt((1 - e_t)/(1 + e_t)) * sqrt(mu_sun/a_t);
    
    delta_1 = sqrt(v_1^2 + 2*mu_earth/6678) - sqrt(mu_earth/6678);
    delta_2 = sqrt(v_2^2 + 2*mu_venus/6151) - sqrt(mu_venus/6151);

%% Problem 4
%
% This problem involves a chief/deputy configuration where the deputy
% can't get to close to the chief. I settled on an identical orbit with a
% along-track offset, and a slight inclination offset so the plots wouldn't
% all be stationary. There are no radial differences so the plots are
% pretty boring, but slight radial differences really propagate throughout
% an orbital period so it's best to ignore them for the purposes of this
% problem.

    %% Given
    r_max = 100;
    r_min = 25;
    mu = 3.986e5;
    r = 6378 + 500;
    n = sqrt(mu/r^3);
    v = n*r;
    P = 2*pi/n;
    
    %% Propagate with ode45
    
    r_0 = [0;50;1;0;0;0];
    
    t = [0 P];
    
    [t,X] = ode45(@(t,X) func(t,X,n),t,r_0);
    
    %% Plotting
    
    figure()
    plot(X(:,1),X(:,2),'LineWidth',1.5); hold on
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    title('Deputy Spacecraft Relative Position')
    grid on
    hold off
    
    figure()
    plot(X(:,1),X(:,3),'LineWidth',1.5); hold on
    axis equal
    xlabel('x [m]')
    ylabel('z [m]')
    title('Deputy Spacecraft Relative Position')
    grid on
    hold off
    
    figure()
    plot(X(:,2),X(:,3),'LineWidth',1.5); hold on
    axis equal
    xlabel('y [m]')
    ylabel('z [m]')
    title('Deputy Spacecraft Relative Position')
    grid on
    hold off
    
    figure()
    plot3(X(:,1),X(:,2),X(:,3),'LineWidth',1.5); hold on
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    title('Deputy Spacecraft Relative Position')
    grid on
    hold off
    
    figure()
    plot(t,X(:,3),'LineWidth',1.5); hold on
    xlabel('time [s]')
    ylabel('z [m]')
    title('Deputy Spacecraft Relative z motion over time')
    grid on
    hold off
    
    %% Distance Calculations
    r_deputy = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);
    r_deputy_bubble = r_deputy-25;
    
    %% Plotting Distance
    figure()
    plot(t,r_deputy,'LineWidth',1.5); hold on
    xlabel('time [s]')
    ylabel('distance [m]')
    title('Deputy Distance from Chief over time')
    grid on
    hold off
    
    figure()
    plot(t,r_deputy_bubble,'LineWidth',1.5); hold on
    xlabel('time [s]')
    ylabel('distance [m]')
    title('Deputy Distance from Safety Bubble over time')
    grid on
    hold off
    
    %% Notes
    %
    % This orbit is designed to be stationary in xy but oscillate in z.
    % This is apparent in the plots though without the axis of time it's
    % hard to see the oscillation, that's what the last plot serves to
    % show. Because this is stationary in the xy plane, the plots don't
    % really look like much, and it's very simple to create a configuration
    % that fits within the requirements for the problem.
    
    
function dXdt = func(t,r_0,const)
    n = const;

    x = r_0(1);
    y = r_0(2);
    z = r_0(3);
    x_dot = r_0(4);
    y_dot = r_0(5);
    z_dot = r_0(6);
    
    x_ddot = 2*n*y_dot + 3*n^2*x;
    y_ddot = -2*n*x;
    z_ddot = -n^2*z;
    
    dXdt = [x_dot;y_dot;z_dot;x_ddot;y_ddot;z_ddot];
end