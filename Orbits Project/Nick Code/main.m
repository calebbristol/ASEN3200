%% House keeping
clc
clear all
close all

%% Variables
% Appearence
E_col = [.2 .2 1]; % Earth's color [0 1 1]
C_col = [.3 1 .3]; % Cities' color [1 0 0]
Co_col = 'k'; % Coastline's color
numCities = 41001;% Number of cities to plot || MAX 41001 ||

% Constants
filename = 'basic.json';
J2 = 1.08263e-3;
mu = 398600.4418; %[km3 s-2]
Re = 6378; %[km]

%% Part 1 -- read in a JSON constellation design file
% [num_launches, num_spacecraft, satellite_list] = loadConstellation(filename);

% satellite_list(1).oe0 = [7103 0.00001 deg2rad(70) deg2rad(110) 0 0];
% satellite_list(2).oe0 = [7103 0.00001 deg2rad(70) deg2rad(105) 0 0];

num_spacecraft = 150;
for i=1:num_spacecraft
    % [ a e i O o f
    satellite_list(i).oe0 = [7103 0.00001 deg2rad(70) i*deg2rad(12/5) 0 0];
end

%% Part 2 -- propagate the constellation in time for a full mean solar day
tstep = 30; %[sec]
for i=1:num_spacecraft
    for j=0:tstep:86400
        unit = j/tstep;
        x(unit+1,:) = propagateState(satellite_list(i).oe0,j,0,mu,J2,Re)';
    end
    satellite_list(i).orbits = x; %Stores orbit position and velocity into satellite_list
    satellite_list(i).r_sc = [x(:,1),x(:,2),x(:,3)]; 
end

%% Part 4 -- plot a 3D rendering of your constellation orbits and the Earth
hold on
axis equal

% Plot Satellite Orbits
for i=1:num_spacecraft
    plot3(satellite_list(i).orbits(:,1),satellite_list(i).orbits(:,2),satellite_list(i).orbits(:,3));
end

% Plot Planet
S = oblateSpheroid; %3D object
S.SemimajorAxis = Re+21; %Make bulge around equation
S.SemiminorAxis = Re;
[xP,yP,zP] = ellipsoid(0,0,0,S.SemimajorAxis,S.SemimajorAxis,S.SemiminorAxis);
h = surf(xP,yP,zP);
set(h,'edgecolor',E_col+[.01 -.1 -.1],'facecolor',E_col)

% Plot Cities
cities = readmatrix('data/worldcities.xlsx');
citylat = cities(1:numCities,3);
citylon = cities(1:numCities,4);
cityPop = cities(1:numCities,10);
[xC,yC,zC] = geodetic2ecef(S, citylat, citylon, 0,'degrees');
scatter3(xC,yC,zC,4,C_col,'.');

% Plot Coastlines
load coastlines
[xL,yL,zL] = geodetic2ecef(S, coastlat, coastlon, 0,'degrees');
plot3(xL,yL,zL,'color',Co_col,'LineWidth',2);
hold off

%% Part 3 -- compute the number of spacecraft in line of sight of each city
r_site=[xC,yC,zC];
elevation_limit = deg2rad(15);

cityAll = 41001;
cityTest = 40;
CITY = cityTest;

timeAll = 2880;
timeTest = 2880;
TIME=timeTest;

info = zeros(TIME,CITY);
% for t=1:TIME
%     for i=1:CITY
%         for j=1:num_spacecraft
%             info(t,i) = info(t,i) + testLoS(r_site(i,:),satellite_list(j).r_sc(t,:),elevation_limit);
%         end
%     end
% end

for j=1:num_spacecraft
    info = zeros(TIME,CITY);
    for t=1:TIME
        for i=1:CITY
            info(t,i) = info(t,i) + testLoS(r_site(i,:),satellite_list(j).r_sc(t,:),elevation_limit);
            
        end
    end
    s(j).info = info;
end
% for t=1:TIME
%     for i=1:CITY
%         j=2;
%         info2(t,i) = info(t,i) + testLoS(r_site(i,:),satellite_list(j).r_sc(t,:),elevation_limit);
%     end
% end

tF = TIME;
t0 = 0;
delT = tF-t0;
delTinv = 1/delT;

% for i=1:CITY
%     part1 = cityPop(i)*delTinv;
%     part2 = sum(info(:,i))*delTinv;
%     Jv(i) = part1*part2;
% end
% JvTotal = sum(Jv);
for j=1:num_spacecraft
    Jv = zeros(1,100);
    for i=1:CITY
        part1 = cityPop(i)*delTinv;
        part2 = sum(s(j).info(:,i));
        Jv(i) = part1*part2;
    end
    JvTotal = sum(Jv);
    s(j).Jtotal = JvTotal;
end
% Jtotal = struct2array(s(:,2));
for p=1:num_spacecraft
    JtotalArray(p) = s(p).Jtotal;
end
[Jmax,pos] = max(JtotalArray);

figure()
plot(1:12/5:360,JtotalArray,'LineWidth',2); hold on
title('Jv Values Based on Varied Omega')
xlabel('Omega [degrees]')
ylabel('Jv')
xlim([0 360])
ylim([0.95*min(JtotalArray) 1.1*max(JtotalArray)])
xline(rad2deg(satellite_list(pos).oe0(4)),'--r','LineWidth',2)
set(gca,'FontSize',14)
legend('Jv Curve','Max Jv Value')
grid on
hold off


