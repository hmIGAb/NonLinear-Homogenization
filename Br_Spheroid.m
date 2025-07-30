clear;close all; clc;
eps0=8.854e-12;
% Material A: Spheroid
a=5; % radius
epsA=4*eps0; % linear 
chiA=9.07571e-12*eps0; % nonlinear
fA=0.5; % volume fraction
% Material B: Spheres
epsB=12*eps0; % linear 
chiB=(0)*eps0; % nonlinear
fB=1-fA; % volume fraction
% Initial value
epsBr0=fA*epsA*eye(3) + fB*epsB*eye(3); % linear 
chiBr0=fA*chiA*eye(3) + fB*chiB*eye(3); % nonlinea

% Bruggeman
tol=1e-5;
E=100;
n=0;
while E>tol
    % disp(n);
    n=n+1;
    % Depolarization dyadics
    [DA0,DA1]=DNL_UX_Spheroid(a,1,epsBr0,chiBr0); % DEP A particles 
    [DB0,DB1]=DNL_UX_Sphere(epsBr0,chiBr0); % DEP B particle
    % Bruggeman function
    [epsBr1,chiBr1]=BrNL_UX(fA,epsA,chiA,DA0,DA1,fB,epsB,chiB,DB0,DB1,epsBr0,chiBr0);
    % Error
    EL =max(abs((diag(epsBr1)-diag(epsBr0))./(diag(epsBr1))));
    ENL=max(abs((diag(chiBr1)-diag(chiBr0))./(diag(chiBr1))));
    E=max([EL,ENL]);
    epsBr0=epsBr1;
    chiBr0=chiBr1;
end
disp(epsBr0/eps0);
disp(chiBr0/chiA);