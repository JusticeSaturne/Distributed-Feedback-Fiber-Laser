function [N1,N2,N3,N5,N6,Gamma,alpha] = rateEquationPairedSolver(Pp,Ps)
% rateEquationsParedSover function solves the Erbium Ytterbium rate equations in
% Steady state taking into account the influence of cooperative
% upconversion the pupulation density of Erbium ions at metastable level
% are divided into two groups (single and paired ions). The algorithm starts with 
% single ions, and paired ions are added by an iterative process
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
%   -Distributed Feedback Fiber Laser Strain and Temperature Sensor, 
%    Olivier Hadeler, University of Southampton 
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
Ners = Ner*0.53;    % Population density of Erbium single ions
Nybs = Nyb*0.53;    % Population density of ytterbium single ions
Nerp = Ner*0.47;    % Population density of Erbium paired ions 
Nybp = Nyb*0.47;    % Population density of Ytterbium paired ions
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
F_s = cel/w; 
den1 = A*hpk*F_s;
% Solve rate equations for single ions
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
g = Ctr*(a*R56*Nybs-b*c*Ners)+e*a*R13+e*b*d;
h = -e*b*c*Ners;
N1s = -(g/(2*f))+((sqrt(g^2-4*f*h))/(2*f));
N3s = (c*Ners-d*N1s)/a;
N2s = Ners-N1s-N3s;
N6s = (R56*Nybs+R13*N1s-b*N3s)/e;
N5s = Nybs-N6s;
Gammasing = Gammas*(sigma21*N2s-sigma12*N1s);
% Gain and loss computation
%----------------solve rate equations for paired ions------------------
icount = 0;
maxit = 1000;
err = 10;
toler = 1e-5;
while (err > toler && icount <= maxit)
    icount = icount+1;
    % Solve rate equations
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
    g = Ctr*(a*R56*Nybp-b*c*Nerp)+e*a*R13+e*b*d;
    h = -e*b*c*Nerp;
    N1p = -(g/(2*f))+((sqrt(g^2-4*f*h))/(2*f));
    N3p = (c*Nerp-d*N1p)/a;
    N2p = Nerp-N1p-N3p;
    N6p = (R56*Nybp+R13*N1p-b*N3p)/e;
    N5p = Nybp-N6p;
    A2new = 100 + Cup*N2p;
    err = abs((A2new-A2)/A2new);
    A2 = A2new;
end
Gammapair = Gammas*(sigma21*N2p-sigma12*N1p);
%------------------Total population densities--------------------------
N1 = N1s+N1p;
N2 = N2s+N2p;
N3 = N3s+N3p;
N5 = N5s+N5p;
N6 = N6s+N6p;
Gamma = Gammas*(sigma21*N2-sigma12*N1)-alphas;
alpha = Gammap*(sigma65*N6-sigma56*N5+sigma31*N3-sigma13*N1)-alphap;