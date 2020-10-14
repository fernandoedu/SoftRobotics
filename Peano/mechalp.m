%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: mechalp.m
% Description:  MATLAB script used for simulating MECHALP
%               static model of a Peano muscle
% Author: Fernando S. Edu Jr.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% Test pressure used in simulation 
Ptest = [0 80000 500000];

% MECHALP model geometric parameters
wt = 0.01725;   % deflated width
lt = 0.046;     % deflated length
alpha = 0.809;  % tube length utilization
hc = 0.002;     % contraction displacement
Ac = hc * lt;   % piston contractile area
la = (pi - 2) * wt / (2 * pi);  % linkage length

e1 = linspace(-0.07, 0, 100);   % contraction ratio range used for 0kPa test
e2 = linspace(0, 0.17, 200);    % contraction ratio range used for 80kPa and 500 kPa tests
Fs = zeros(1, size(e1,2) + 2 * size(e2,2)); % Force solved for MECHALP model
Fc = zeros(1, 2 * size(e2,2));   % Force solved for Niiyama model

Ce = 1.0934*10^-6;     % Niiyama's correction factor

off_1 = size(e1,2);                 % start offset for 80kPa data
off_2 = size(e1,2) + size(e2,2);    % start offset for 500kPa data

% select levenberg-marquardt algorithm for optimization
opt = optimoptions('fsolve');
opt.Algorithm = 'levenberg-marquardt';

for i = 1:size(Fs,2)
    if(i <= off_1)
        % simulate 0 kPa set-up for MECHALP model only
        P = Ptest(1);
        e = e1;
        
        xsol = fsolve(@(x) P * (0.91 * alpha * wt * lt * (la - ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + (e(i)*wt)) / 2)) ...
            - sqrt(la^2 - (la - ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + (e(i)*wt)) / 2)^2) ...
            * ((P * Ac) + (((1*10^15) * ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + e(i)*wt)^4) * e(i) * wt) ...
            + (1 + (((1*10^15) * ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + e(i)*wt)^4) / (0.5 ...
            * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x)))))) * x), 500, opt);
        Fs(i) = xsol;
        
    elseif((i > off_1) && (i <= off_2))
        % simulate 80 kPa set-up for both MECHALP and Niiyama models
        P = Ptest(2);
        e = e2;
        
        % solve force for MECHALP model
        xsol = fsolve(@(x) P * (0.91 * alpha * wt * lt * (la - ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + (e(i - off_1)*wt)) / 2)) ...
            - sqrt(la^2 - (la - ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + (e(i - off_1)*wt)) / 2)^2) ...
            * ((P * Ac) + (((1*10^15) * ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + e(i - off_1)*wt)^4) * e(i - off_1) * wt) ...
            + (1 + (((1*10^15) * ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + e(i - off_1)*wt)^4) / (0.5 ...
            * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x)))))) * x), 500, opt);
        Fs(i) = xsol;
        
        % solve force for Niiyama model
        tsol = fsolve(@(x) (1 - (sin(x) / x)) * (1 + ((pi / (pi - 2)) * Ce * P)) - (Ce * P) - e(i - off_1), 1);
        Fc(i - off_1) = P * wt * lt * (cos(tsol) / tsol);
    else
        % simulate 500 kPa set-up for both MECHALP and Niiyama models
        P = Ptest(3);
        e = e2;
        
        % solve force for MECHALP model
        xsol = fsolve(@(x) P * (0.91 * alpha * wt * lt * (la - ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + (e(i - off_2)*wt)) / 2)) ...
            - sqrt(la^2 - (la - ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + (e(i - off_2)*wt)) / 2)^2) ...
            * ((P * Ac) + (((1*10^15) * ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + e(i - off_2)*wt)^4) * e(i - off_2) * wt) ...
            + (1 + (((1*10^15) * ((x/(0.5 * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x))))) + e(i - off_2)*wt)^4) / (0.5 ...
            * (8*10^4 + sqrt(6.4*10^9 + (6*10^8 * x)))))) * x), 500, opt);
        Fs(i) = xsol;
        
        % solve force for Niiyama model
        tsol = fsolve(@(x) (1 - (sin(x) / x)) * (1 + ((pi / (pi - 2)) * Ce * P)) - (Ce * P) - e(i - off_2), 1);
        Fc(i - off_1) = P * wt * lt * (cos(tsol) / tsol);
    end
end

% Create object for figure and axes
fh = figure;
ah = axes(fh);

hold on;

% Plot contractile strain vs. force
plot(ah,e1*100,Fs(1:100),'bo','DisplayName','MECHALP @ 0kPa');
plot(ah,e2*100,Fs(101:300),'kx','DisplayName','MECHALP @ 80kPa');
plot(ah,e2*100,Fc(1:200),'r--','DisplayName','Niiyama @ 80kPa');
plot(ah,e2*100,Fs(301:500),'gs','DisplayName','MECHALP @ 500kPa');
plot(ah,e2*100,Fc(201:400),'md','DisplayName','Niiyama @ 500kPa');
xlabel('Contractile strain (%)');
ylabel('Force (N)');

axis(ah,[-7,20,0,350])
grid on;
legend show;
