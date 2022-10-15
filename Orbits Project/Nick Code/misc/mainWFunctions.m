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
filename = 'example_constellation.json';
J2 = 1.08263e-3;
mu = 398600.4418; %[km3 s-2]
Re = 6378; %[km]
tstep = 30; %[sec]

%% Part 1 -- read in a JSON constellation design file
[num_launches, num_spacecraft, satellite_list] = loadConstellation(filename);

%% Part 2 -- propagate the constellation in time for a full mean solar day
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
cities = readmatrix('worldcities.xlsx');
citylat = cities(1:numCities,3);
citylon = cities(1:numCities,4);
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
info = zeros(2880,41001);
for i=1:2880
    for j=1:41%001
        for k=1:num_spacecraft
            info(i,j) = info(i,j) + testLoS(r_site(j,:),satellite_list(k).r_sc(i,:),elevation_limit);
        end
    end
end

%% Functions

function [num_launches, num_spacecraft, satellite_list] = loadConstellation(filename)
    num_launches = 0;
    num_spacecraft = 0;
    satellite_list.name = '';
    satellite_list.oe0 = NaN(6,1);

    %1) extract the constellation structure from the json file
    fid = fopen(filename);
    raw = fread(fid,inf);
    str = char(raw'); 
    fclose(fid); 
    data = jsondecode(str);

    %2) read all of the launches and payloads to understand how many launches
    % and spacecraft are in the constellation; note, this will be useful in
    % Part 2!

    num_launches = length(data.launches);
    for i=1:num_launches
        num_spacecraft = num_spacecraft + length(data.launches(i).payload);
    end

    %3) RECOMMENDED: Pre-allocate the satellite_list struct
    satellite_list(num_spacecraft).name = '';
    satellite_list(num_spacecraft).oe0 = NaN(6,1);

    %4) Populate each entry in the satellite struct list with its name and
    %initial orbit elements [a,e,i,Om,om,f] at time t0
    i = 1;
    for j=1:num_launches
        oe0Temp1 = cell2mat(struct2cell((data.launches(j).orbit)));
        for k=1:length(data.launches(j).payload)
            nameTemp = data.launches(j).payload(k).name;
            oe0Temp2 = data.launches(j).payload(k).f;
            satellite_list(i).name = nameTemp;
            satellite_list(i).oe0 = [oe0Temp1; oe0Temp2];
            i = i + 1;
        end
    end
end

function x = propagateState(oe0,t,t_0,mu,J2,Re)
    %make sure that function has outputs
    x = NaN(6,1);

    %misc;
    a = oe0(1);
    e = oe0(2);
    i = oe0(3);
    Om = oe0(4);
    om = oe0(5);
    f = oe0(6);

    p = a*(1-e^2);
    h = sqrt(mu*p);

    %1) Compute the mean orbit elements oe(t) at time t due to J2 perturbations
    n = sqrt(mu/a^3); %Mean motion
    M = n*(t-t_0); %Mean anomoly

    OmDot = (-3/2)*n*J2*(Re/p)^2*cos(i); %OMEGA dot
    omDot = (3/2)*n*J2*(Re/p)^2*(2-(5/2)*sin(i)^2); %omega dot
    dOm = (t-t_0)*OmDot; %Change (delta) in OMEGA
    dom = (t-t_0)*omDot; %...omega
    Om = Om + dOm; %New OMEGA from J2 perturbations
    om = om +dom; %...omega

    %2) Solve the time-of-flight problem to compute the true anomaly at time t

    f = M+(2*e-.25*e^3)*sin(M)+(5/4)*e^2*sin(2*M)+(13/12)*e^3*sin(3*M);

    %3) Compute r(t), rdot(t) in the perifocal frame
    r = p/(1+e*cos(f));
    rVecP = r*[cos(f); sin(f); 0];
    vVecP = (mu/h)*[-sin(f); e+cos(f); 0];

    %4) Compute r(t), rdot(t) in the ECI frame, save into x
    DCM_EP = toEfromP(Om,i,om);
    rVecE = DCM_EP * rVecP;
    vVecE = DCM_EP * vVecP;
    x = [rVecE; vVecE];
end

function inLoS = testLoS(r_site,r_sc,elevation_limit)
    foundAngle = dot(r_site,r_sc)/(norm(r_site)*norm(r_sc));

    if foundAngle > elevation_limit
        inLoS = 1;
    else
        inLoS = 0;
    end
end

function DCM = toEfromP(OMEGA,i,omega)
    cO = cos(OMEGA);
    sO = sin(OMEGA);
    ci = cos(i);
    si = sin(i);
    co = cos(omega);
    so = sin(omega);
    DCM = [co*cO-so*ci*sO co*sO+so*ci*cO so*si;...
            -so*cO-co*ci*sO -so*sO+co*ci*cO co*si;...
            si*sO -si*cO ci]';
end
