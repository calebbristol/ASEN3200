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