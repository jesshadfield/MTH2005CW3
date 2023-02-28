close all
clear  
clc

% Cosntants as given in background info  
g = 9.81;                        % Acceleration due to gravity (m s^-2)
c_pa = 1005.0;                   % Specific heat capacity of dry air (J kg^-1 K^-1)
Rho_w = 1000.0;                  % Density of liquid water (Kg m^-3)
Rho_a = 1.225; 			         % Density of air (Kg m^-3)
Eps = 0.622;                     % Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6;                      % Latent heat of vapourisation (J Kg^-1)
Ra = 287.0;                      % Gas constant of dry air (J kg^-1 K^-1)
Rv = 462.0;                      % Gas constant of water vapour (J kg^-1 K^-1)
k = 0.024;                       % Thermal Conductivity of Air (J m^-1 s^-1 K^-1)
Kv = 2.21e-5;                    % Diffusivity of Water Vapour (m^2 s^-1)


X = [0.0035, 280]; % X = [Supersaturation, Temp(K)] 

[es] = svp(X);	%! Saturation vapor pressure (provided by function svp.m)


% Formula for A3
A3 = (((Lv)^2 * Rho_w)/(k*Rv*(X(2))^2) + (Rho_w * Rv * X(2))/(Kv*es))^(-1);

s = X(1); % Supersaturation 

% Initial Radius in meters squared 
r = (0.85 * 1e-6); 

% Time step set up

%deltat = (1/60); 
deltat = 1;

t = 0:deltat:45; % 0 to 45 minutes broken up into sections 

nstop = length(t) - 1; % Setting up nstop 

% Forward Euler 

for n = 1:nstop 
r(n+1) = r(n) + BLCDG(r(n), s, A3)*deltat;
end 

% Plot

figure(1)
plot (t*1e+6,r, 'r')
xlabel('Time')
ylabel('Radius of Droplet')
set(gca,'fontsize',14)


% Fourth Order Runge-Kutta 

k1 = BLCDG(r(n), s, A3);
k2 = BLCDG(r(n)+0.5*deltat*k1, s, A3);
k3 = BLCDG(r(n)+0.5*deltat*k2, s, A3);
k4 = BLCDG(r(n)+deltat*k3,s, A3);
r(n+1) = r(n) + (deltat/6)*(k1+2*k2+2*k3+k4);

figure(2)
plot(t*1e+6, r, 'b')
xlabel('Time')
ylabel('Radius of Droplet')
set(gca,'fontsize',14)



