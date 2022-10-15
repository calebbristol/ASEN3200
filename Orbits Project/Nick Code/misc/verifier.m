%DEVIN WU
clear
clc

%The purpose of this code is to assist users with aeronautic calculations.

% acquire user inputs
disp('Format of entering vectors: [x,y,z]');
position=input('Enter the position vector in ECI coordinates:');
velocity=input('Enter the velocity vector in ECI coordinates:');
% t1=input('Enter an initial time (sec):');
% t2=input('Enter a later time (sec):');
t1 = 0;
t2 = 10;
% relative=input('Enter a relative tolerance value: ');
% absolute=input('Enter an  absolute tolerance value: ');
relative = 1e-9;
absolute = relative
% mu=398600; %KM
mu = 1;
I=[1,0,0];
J=[0,1,0];
K=[0,0,1];

%output the input data
disp(' ');
disp('Position vector:');
disp(position);
disp('Velocity vector: ');
disp(velocity);
% disp('Initial time (sec): ');
% disp(t1);
% disp('Later time (sec): ');
% disp(t2);
% disp('Relative tolerance: ');
% disp(relative);
% disp('Absolute tolerance: ');
% disp(absolute);

%calculation
r=norm(position);
v=norm(velocity);
energy=(v^2)/2-mu/r;
a=-mu/(2*energy);

angularmo=cross(position,velocity);
h=norm(angularmo);
evector=cross(velocity,angularmo)/mu-position/r;
e=norm(evector);

i=acos(dot(angularmo,K)/h);

n=cross(K,angularmo)/norm(cross(K,angularmo));
if dot(I,n)>0
    omega=atan(dot(J,n)/dot(I,n));
else
    omega=atan(dot(J,n)/dot(I,n))+pi;
end

if dot(evector,K)>0
    w=acos(dot(n,evector)/e);
else
    w=-acos(dot(n,evector)/e);
end

c3w=[cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1];
c1i=[1,0,0;0,cos(i),sin(i);0,-sin(i),cos(i)];
c3omega=[cos(omega),sin(omega),0;-sin(omega),cos(omega),0;0,0,1];
cep=c3w*c1i*c3omega;

transrt1=cep*transpose(position);
rt1=transpose(transrt1);
transvt1=cep*transpose(velocity);
vt1=transpose(transvt1);

p=h^2/mu;
T=2*pi*sqrt(a^3/mu);
if vt1(1)<0
    theta1=acos((p/r-1)/e);
else 
    theta1=-acos((p/r-1)/e);
end
Et1=2*atan(sqrt((1-e)/(1+e))*tan(theta1/2));

Mt1=Et1-e*sin(Et1);
T0=t1-Mt1/sqrt(mu/a^3);

Mt2=sqrt(mu/a^3)*(t2-T0);
Eold=1;
Enew=Mt2;
while or(abs(Enew-Eold)>absolute,abs((Enew-Eold)/Eold)>relative)
    Eold=Enew;
    Enew=Eold-(Eold-e*sin(Eold)-Mt2)/(1-e*cos(Eold));
   
end
Et2=Enew;
thetat2=2*atan(sqrt((1+e)/(1-e))*tan(Et2/2));

r2=a*(1-e*cos(Et2));
rp2=[r2*cos(thetat2),r2*sin(thetat2),0];
vp2=[-sqrt(mu/p)*sin(thetat2),sqrt(mu/p)*(e+cos(thetat2)),0];

transrp2=transpose(cep)*transpose(rp2);
position2=transpose(transrp2);

transvp2=transpose(cep)*transpose(vp2);
velocity2=transpose(transvp2);

%output computation results
disp('a :');
disp(a);
disp('e :');
disp(e);
disp('i :');
disp(i);
disp('OMEGA :');
disp(omega);
disp('omega :');
disp(w);
disp('f t1 :');
disp(theta1);
% disp('DCM cep :');
% disp(cep);
disp('The position vector (in perifocal coordinates) at time t1 (km):');
disp(rt1);
disp('The velocity vector (in perifocal coordinates) at time t1 (km/s):');
disp(vt1);
% disp('The eccentric anomaly at time t1 (rad):');
% disp(Et1);
% disp('The time of periapsis passage (sec):');
% disp(T0);
% disp('The true anomaly at time t2 (rad):');
% disp(thetat2);
% disp('The eccentric anomaly at time t2 (rad):');
% disp(Et2);
% disp('The position vector (in perifocal coordinates) at time t2 (km):');
% disp(rp2);
% disp('The velocity vector (in perifocal coordinates) at time t2 (km/s):');
% disp(vp2);
% disp('The position vector (in ECI coordinates) at time t2 (km):');
% disp(position2);
% disp('The velocity vector (in ECI coordinates) at time t2 (km/s):');
% disp(velocity2);