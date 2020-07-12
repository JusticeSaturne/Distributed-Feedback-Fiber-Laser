function [N1,N2,N3,N5,N6,Gamma,alpha] = rateEquationSolver(Pp,Ps)
% rateEquationsSover function solves the Erbium Ytterbium rate equations in
% Steady state without taking into account the influence of cooperative
% upconversion
%
% Inputs:
% - Pp : Pump power
% - Ps: Laser power
%
% Outputs:
%   N1: Population density of the ground state level of Erbium
%   N2: Population density of the metastable level of Erbium
%   N3: Population density of the energy level 3 of Erbium
%   N5: Population density of the ground state of Ytterbium 
%   N6: Population density of the upper level of Ytterbium
%   Gamma: Gain coefficient of the laser
%   alpha: Absorption coefficient of the propagating field inside the
%   cavity
%
% Comments:
%   - Computation is done in the ideal case of a two level energy system
%
% References:
%
% Written by Justice Sompo, University of Johannesburg, South Africa
%
% Rate equation Parameters
sigma12 = 0.28e-24;     % Absoption cross section of Erbium ions at 1550 nm (m^2)
sigma21 = 0.42e-24;     % Emission cross section of Erbium ions at 1550 nm (m^2)
sigma56 = 0.58e-24;     % Absorption cross section of Ytterbium ions at 980 nm (m^2)    
sigma65 = 2.2e-24;      % Emission cross section of Ytterbium ions at 980 nm (m^2) 
sigma13 = 0.20e-24;     % Absorption cross section of Erbium at 980 nm (m^2)
sigma31 = 0.20e-24;     % Emission cross section of Erbium at 980 nm (m^2)
% Lifetimes
t65 = 0.8e-3;           % Lifetime of upper energy level of Ytterbium
t21 = 10e-3;            % Lifetime of metastable energy level of Erbium
t31 = 10e-6;            % Lifetime of upper energy level of Erbium
A6 = 1/t65;
A3 = 1/t31;
A2 = 1/t21;
% ion-ions interactions
Cup = 2.5e-21;      % cooperative upconversion coefficient
Ctr = 5e-21;        % Er-Yb energy transfer coefficient
% Concentrations
Ner = 1.2e25;       % Population density of Erbium ions
Nyb = 24e25;        % Population density of Ytterbium ions
% Overlap factors and background losses
Gammap = 0.82;      % Erbium overlap factor at pump wavelength
Gammas = 0.73;      % Erbium overlap factor at laser wavelength
alphap = 0.2;       % Background loss of the Er-Yb fibre at 980 nm
alphas = 0.15;      % Background loss of the Er-Yb fibre at laser wavelength of 1550 nm
% pump wavelength
lambdap = 980e-9;
w =1547.3661e-9; % Output wavelength
% Other physical constants
r = 2.3e-6;
hpk = 6.626e-34;
cel = 3e8;
% Calculated parameters
A = pi*r^2;
F_p = cel/lambdap;
den2 = A*hpk*F_p;
F_s = cel/w;den1 = A*hpk*F_s;
W12 = (sigma12*Ps*Gammas)/den1;
W21 = (sigma21*Ps*Gammas)/den1;
R56 = (sigma56*Pp*Gammap)/den2;
R65 = (sigma65*Pp*Gammap)/den2;
R13 = (sigma13*Pp*Gammap)/den2;
a = W21+A2+A3;
b = R13+A3;
c = W21+A2;
d = W12+W21+A2;
e = R56+R65+A6;
f = Ctr*(a*R13+b*d);
g = Ctr*(a*R56*Nyb-b*c*Ner)+e*a*R13+e*b*d;
h = -e*b*c*Ner;
N1 = -(g/(2*f))+((sqrt(g^2-4*f*h))/(2*f));
N3 = (c*Ner-d*N1)/a;
N2 = Ner-N1-N3;
N6 = (R56*Nyb+R13*N1-b*N3)/e;
N5 = Nyb-N6;
% Gain and loss computation
Gamma = Gammas*(sigma21*N2-sigma12*N1)-alphas;
alpha = Gammap*(sigma65*N6-sigma56*N5+sigma31*N3-sigma13*N1)-alphap;