%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: niiyama.m
% Description:  MATLAB script used for simulating Niiyama's
%               model of a Pouch motors
% Author: Fernando S. Edu Jr.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

P = [0.25 0.5 1]';  % Test pressure used in simulation 

% Peano muscle geometric properties
wt = 1;
lt = 1;

theta = [0.01:0.01:pi/2];   % range of included angle os muscle sector

% Ideal Peano model
ei = (1 - (sin(theta) ./ theta)) .* 100;
Fi = P * wt * lt * (cos(theta) ./ theta);

% Niiyama et. al's Pouch motor model
Ce = 0.025;     % corrective coefficient
ec = ((1 - (sin(theta) ./ theta)) .* (1 + ((pi / (pi - 2)) * Ce * P)) - (Ce * P)) .* 100;
Fc = P * wt * lt * (cos(theta) ./ theta);

% Create object for figure and axes
fh = figure;
ah = axes(fh);

hold on;

% Plot contractile strain vs. force
plot(ah,ei,Fi(1,:),'ro','DisplayName','ideal P = 0.25');
plot(ah,ei,Fi(2,:),'g+','DisplayName','ideal P = 0.5');
plot(ah,ei,Fi(3,:),'b*','DisplayName','ideal P = 1');
plot(ah,ec(1,:),Fc(1,:),'c.','DisplayName','corrected P = 0.25');
plot(ah,ec(2,:),Fc(2,:),'mx','DisplayName','corrected P = 0.5');
plot(ah,ec(3,:),Fc(3,:),'ks','DisplayName','corrected P = 1');
xlabel('Contractile strain (%)');
ylabel('Force (N)');

axis(ah,[0,35,0,4])
grid on;
legend show;

