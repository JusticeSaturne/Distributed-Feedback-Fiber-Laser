function [S,GAMMA,PP,PS,N1,N2,N3,N5,N6] = dfbSimulation(Pp,P0)
% Function dfbsimulation compute the DFB fibre laser internal gain, pump
% power and laser power distribution given pump power and laser power
% 
% 
% Inputs:
%  - Pp:Pump power 
%  - P0:Laser power 
%  
% Outputs: 
%   - N1: Population density distribution at ground level of Erbium
%   - N2: Population density distribution at metastable level of Erbium
%   - N3: Population density distribution at upper level of Erbium
%   - N5: Population density distribution at ground level of Ytterbium
%   - N6: Population density distribution at upper level of Erbium
%   - Gain distribution inside the cavity
%   - PP: Pump power distribution
%   - PS: Laser power distribution inside the cavity
%   - GAMMA: Gain distribution
% Comments:
%   - Erbium-Ytterbium co-doped gain medium, rate equation solved
%   analytically neglecting cooperative upconversion
%
% References: 
%  
% Written by Justice Sompo, University of Johannesburg, South Africa   
format longe
% Fiber Bragg Grating Parameters
L = 0.050;M = 100;
dn =7.38e-5;
neff= 1.47; 
lambdaB = 1547.3661e-9; 
period = lambdaB/(2*neff);
dz = L/M;
phaseshift =[exp(-1i*(pi/2)) 0;0 exp(1i*(pi/2))];
kac = (pi*dn)/lambdaB; 
phase =0;
% Initial guess
w =1547.3661e-9;
%Pp = 100e-3;
%P0 =18.13e-3;
S0 = sqrt(P0); %Estimated backward laser amplitude
R0 = 0;       % No forwards propagating wave
F = [R0;S0];
Ps = (abs(S0))^2+(abs(R0))^2;
N1 = zeros(1,100);
for ii=1:M
    [n1,n2,n3,n5,n6,Gamma,alpha] = rateEquationSolver(Pp,Ps);
    % finish solving the rate equation here
    %Pump power distribution
    Pp = Pp*exp(alpha*dz);
    % Feedback mechanism
    beta = pi/period;
    dets = 2*pi*neff*(1/w-1/lambdaB);
    detsprime = (dets+1i*Gamma);
    p1 = sqrt(kac^2 - detsprime^2);
    % Transfer matrix
    f11 = (cosh(p1*dz)-1i*((dets+1i*Gamma)/p1)*sinh(p1*dz))*exp(-1i*beta*dz);
    f12 = (kac/p1)*sinh(p1*dz)*exp(-1i*(beta*dz+phase));
    f21 = (kac/p1)*sinh(p1*dz)*exp(1i*(beta*dz+phase));
    f22 = (cosh(p1*dz)+1i*((dets+1i*Gamma)/p1)*sinh(p1*dz))*exp(1i*beta*dz);
    ff = [f11 f12;f21 f22];
    F = ff*F;
    if ii == 50
        F = phaseshift*F;
    end
    phase = phase+(2*pi*dz)/period; 
    R=F(1,1);
    S=F(2,1);
    Ps = (abs(S))^2+(abs(R))^2;
    PLb = (abs(S))^2;
    PLf = (abs(R))^2;
    N1(ii)=n1;
    N2(ii)=n2;
    N3(ii)=n3;
    N5(ii)=n5;
    N6(ii)=n6;
    GAMMA(ii)=Gamma;
    PP(ii)=Pp;
    PS(ii)=Ps;
end
