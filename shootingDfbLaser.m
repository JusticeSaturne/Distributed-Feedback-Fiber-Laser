function [S,Ps] = shootingDfbLaser(Pp)
% ShootingDfbLaser function computes the output power of a DFB fibre laser
% using a shooting algorithm with two initial guesses and a secant method
% for guess refinment
% 
% Inputs:
% - Pp : Pump power
% 
%  
% Outputs: 
%   S: Forward internal eletric field
%   Ps: Forward output power  
% Comments:
%   - The function calls the function dfbSimulation
%   - The function " dfbSimulation" is called multiple time to compute
%   different values like:
%   - N1: Population density distribution at ground level of Erbium
%   - N2: Population density distribution at metastable level of Erbium
%   - N3: Population density distribution at upper level of Erbium
%   - N5: Population density distribution at ground level of Ytterbium
%   - N6: Population density distribution at upper level of Erbium
%   - Gain distribution inside the cavity
%   - PP: Pump power distribution
%   - PS: Laser power distribution inside the cavity
%
% References: 
%  
% Written by Justice Sompo, University of Johannesburg, South Africa 
% 
%--------------Initial guess to initiate the Shooting algorithm------------
P01 = 1e-3;
P02 = 10e-3;
L = 50; % length of the cavity
ii = 1;
S1 = dfbSimulation(Pp,P01); % final backward propagation wave for initial guess P01
S2 = dfbSimulation(Pp,P02); % final backward propagation wave for initial guess P02
while(abs(P01-P02)>0.000001)
    tmp = P02;
    P02 = P01-(P01-P02)/(S1-S2)*S1;
    P01 = tmp;
    S1 = dfbSimulation(Pp,P01); % final backward propagation wave for initial guess P01
    S2 = dfbSimulation(Pp,P02); % final backward propagation wave for initial guess P02
    ii = ii+1;
end
ii;
S = S2;
Ps = abs(P02);
%-----------------------------------PLOTTING-------------------------------
[S,GAMMA,PP,PS,N1,N2,N3,N5,N6] = feval('dfbSimulation',Pp,Ps);
z = linspace(0,L);
figure(1)
subplot(2,2,1)
plot(z,N2,'r','Linewidth',2)
xlabel('Position (mm)')
ylabel('N2 Population (ions/m^3)')
subplot(2,2,2)
plot(z,GAMMA,'m','Linewidth',2)
xlabel('Position (mm)')
ylabel('Gain (m^-1)')
subplot(2,2,3)
plot(z,PP,'g','Linewidth',2)
xlabel('Position (mm)')
ylabel('Pump Power (Watts)')
subplot(2,2,4)
plot(z,PS,'Linewidth',2)
xlabel('Position (mm)')
ylabel('Laser Power (Watts)')
figure(2)
subplot(2,2,1)
plot(z,N1,'r','Linewidth',2)
xlabel('Position (mm)')
ylabel('N1 Population (ions/m^3)')
subplot(2,2,2)
plot(z,N3,'m','Linewidth',2)
xlabel('Position (mm)')
ylabel('N3 Population (ions/m^3)')
subplot(2,2,3)
plot(z,N5,'g','Linewidth',2)
xlabel('N5 Population (ions/m^3)')
ylabel('Pump Power (Watts)')
subplot(2,2,4)
plot(z,N6,'Linewidth',2)
xlabel('Position (mm)')
ylabel('N6 Population (ions/m^3)')