%%%%%%%%%%%
% MATLABBIN
%
%%%%%%%%%%%

%% 
clc
clear
close all;

%%
K_p = linspace(0,200,100);
K_d = linspace(0,30,100);

%%
omega_n = @(K_p) sqrt(K_p);
xi = @(K_d,K_p) K_d / (2 * omega_n(K_p));
omega_d = @(xi,omega_n) omega_n * sqrt(1 - xi^2);

%%
Beta = @(xi) atan(sqrt((1 - xi^2) / xi));

t_r = @(xi) (pi - Beta(xi)) / omega_d;
t_p = @(omega_d) pi ./ omega_d;
M_p = @(xi) exp(-(pi / tan(Beta(xi))));
t_s = @(omega_n,xi,p) -log(p/100)/(omega_n * xi);

yes = zeros(100);

for i = 1:100
    for j = 1:100
        o_n = omega_n(K_p(i));
        xi_ = xi(K_d(j),K_p(i));
        o_d = omega_d(xi_,o_n);
        
        %t_r_ = t_r(xi_);
        %t_p_ = t_p(o_d);
        M_p_ = M_p(xi_);
        t_s_ = t_s(o_n,xi_,5);
        
        if t_s_ < 1.5 && M_p_ < 0.1
            yes(i,j) = 1;
        end
    end
end


