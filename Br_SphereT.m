clear;close all; clc;
eps0=8.854e-12;
% Material A: Truncated Spheres
k=0.5;
epsA=4*eps0; % linear 
chiA=9.07571e-12*eps0; % nonlinear
fA=0.5; % volume fraction
% Material B: Spheres
epsB=12*eps0; % linear 
chiB=(0)*eps0; % nonlinear
fB=1-fA; % volume fraction
% Initial value
epsBr0=(fA*epsA + fB*epsB)*eye(3); % linear 
chiBr0=(fA*chiA + fB*chiB)*eye(3); % nonlinear

% Bruggeman
E=100;
n=0;
tol=1e-5;
while E>tol
    n=n+1;
    % Depolarization dyadics
    [DA0,DA1]=DNL_UX_SphereT_INT(k,epsBr0,chiBr0); % DEP A particles 
    [DB0,DB1]=DNL_UX_Sphere(epsBr0,chiBr0); % DEP A particles 
    % Bruggeman function
    [epsBr1,chiBr1]=BrNL_UX(fA,epsA,chiA,DA0,DA1,fB,epsB,chiB,DB0,DB1,epsBr0,chiBr0);
    % Error
    EL =max(abs((diag(epsBr1)-diag(epsBr0))./(diag(epsBr1))));
    ENL=max(abs((diag(chiBr1)-diag(chiBr0))./(diag(chiBr1))));
    E=max([EL,ENL]);
    % disp([n,E]);
    epsBr0=epsBr1;
    chiBr0=chiBr1;
end
disp(epsBr0/eps0);
disp(chiBr0/chiA);