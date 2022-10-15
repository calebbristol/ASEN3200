%% ASEN 3200 Orbital Project Part 1
%  
% Problem Statement:
%
% Author: Caleb Bristol
% Collaborators: N/A
% Date: 12/01/21
%      


%% Workspace Cleaning

clc
clear
close all;

%% Preliminary Work
%
% Stuff not specified in any part but required for calculations.

    %% Constants
    mu = 3.986e14 * 1e-9; %[km^3/s^2] 1e-9 to convert from m^3 to km^3
    Re = 6378.137; %[km]
    J2 = 1082.63*10^(-6);

    %% Read in Cities & Coastline 
    %
    % Used for plotting later need conversion to cartesian for 3D
    cities_data = readtable('worldcities.csv');
    cities_ll = [cities_data.lng cities_data.lat];
    coastlines_ll = load('world_coastline_low.txt');
    
    [cities(:,1),cities(:,2),cities(:,3)] = sph2cart(deg2rad(cities_ll(:,1)),deg2rad(cities_ll(:,2)),Re);
    
    [coastlines(:,1),coastlines(:,2),coastlines(:,3)] = sph2cart(deg2rad(coastlines_ll(:,1)),deg2rad(coastlines_ll(:,2)),Re);

%% Determine constellations
%
% Create a json file with constellations of varying orbital elements

    %% Radius
    %
    % Re + 350 <= Rp <= Ra <= Re + 1100
    
    %% 