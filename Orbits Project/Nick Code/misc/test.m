% for i=1:30:86401
%     i
% end

% [X,Y,Z] = sphere;
% surf(5*X,5*Y,5*Z)
% axis equal

% figure
% for i = 1:40
%     plot3(i,abs(sqrt(100-i^2)),5,'o');
%     xlim([-10 10])
%     ylim([-10 10])
%     zlim([-10 10])
%     hold on
%     pause(.1)
% end 
% hold off
% for i = 1:40
%     [r,v] = orb2rv(10000,1.0,pi/2,0,0,0,0,0,0);
% fdot = pi/4

% tstep=30;
% for i=0:tstep:86400
%     i/tstep
% %     [i*tstep i*tstep-tstep]
% end

% r_site=[6371, 6371, 0]';
% r_sc=[7000, 7000, 0]';
% ele_lim=pi/8;
% yReference=testLoS(r_site, r_sc, ele_lim)

% Load coastline data

% axesm hatano
% meshm(topo60c,topo60cR)
% zlimits = [min(topo60c(:)) max(topo60c(:))];
% demcmap(zlimits)
% colorbar

% % Initialize geoglobe
% f = uifigure();
% g = geoglobe(f);
% % Plot coastlines on geoglobe
% geoplot3(g, coastlat, coastlon, 0, "Color", "red", "LineWidth", 5);]

clear all 
close all
clc 

hold on
Re = 6378; %km
% map = [0 1 1];
% colormap(map)
load('topo60c')
demcmap(topo60c)

S = oblateSpheroid;
S.SemimajorAxis = Re+21;
S.SemiminorAxis = Re;
[xP,yP,zP] = ellipsoid(0,0,0,S.SemimajorAxis,S.SemimajorAxis,S.SemiminorAxis);
h = surf(xP,yP,zP);
set(h,'edgecolor','k','facecolor','w')

load coastlines
[xL,yL,zL] = geodetic2ecef(S, coastlat, coastlon, 0,'degrees');
plot3(xL,yL,zL,'color','k','LineWidth',2);
axis equal

tp = load('topo60c')
% tp = tp + Re;
% % geoshow(topo60c,topo60cR,"DisplayType","texturemap")
% demcmap(topo60c+Re)
% [xL,yL,zL] = geodetic2ecef(S, coastlat, coastlon, 0,'degrees');
% plot3(xL,yL,zL,'color','k','LineWidth',2);
% worldmap(topo60c,topo60cR)
% 
% geoshow(topo60c,topo60cR,'DisplayType','texturemap')
% demcmap(topo60c)
% load topo60c
% axesm hatano
% meshm(topo60c,topo60cR)
% zlimits = [min(topo60c(:)) max(topo60c(:))];
% demcmap(zlimits)
% colorbar