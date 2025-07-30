function [D_L,D_NL]=DNL_UX_SphereT_INT(k,epsL,chi)

% clear; clc;
% epsL = diag([[1,1]*rand*5,rand*5]);
% chi=diag([[1,1]*rand*5,rand*5])*1e-4;
% k=rand;

epstL=epsL(1,1);
epszL=epsL(3,3);
% Components of epsNL
chit=chi(1,1);
chiz=chi(3,3);
% gamma
gammaL=epszL/epstL;
gammaNL=chiz/epstL-(epszL*chit/epstL^2);

A=diag(1./diag(epsL));
B=-A*chi*A;

trOC1=@(th,Ct,Cz) Ct*sin(th).^2 + Cz*(cos(th)-(1-k)).^2;
trOC2=@(rho,Ct,Cz) Ct*rho.^2 + Cz*k^2;

dS1=@(th) sin(th);
% Linear 
LtL=1/(4) * integral(@(th) ... % Transverse component > NonLinear
    sin(th).^2./...
    ((trOC1(th,A(1,1),A(3,3))*epszL).^(3/2)).*dS1(th), ...
    0,acos(1-2*k));
LzL1=1/(2) * integral(@(th) ...% Z component > NonLinear
    (cos(th)).*(cos(th)-(1-k))./...
    ((trOC1(th,A(1,1),A(3,3))*epszL).^(3/2)).*dS1(th), ...
    0,acos(1-2*k));
LzL2=1/(2)*integral(@(rho) ...
    k*rho./...
    (trOC2(rho,A(1,1),A(3,3))*epszL).^(3/2),...
    0,2*sqrt(k*(1-k)));
LzL=LzL1+LzL2;

L_L=diag([LtL,LtL,LzL]);
% NonLinear
LtNL=-3/(8) * integral(@(th) ... % Transverse component > NonLinear
    sin(th).^2.*(chiz*trOC1(th,A(1,1),A(3,3)) + epszL*trOC1(th,B(1,1),B(3,3)))./...
    ((trOC1(th,A(1,1),A(3,3))*epszL).^(5/2)).*dS1(th), ...
    0,acos(1-2*k));
LzNL1=-3/(4) * integral(@(th) ...% Z component > NonLinear
    ((cos(th)).*(cos(th)-(1-k)).*(chiz*trOC1(th,A(1,1),A(3,3)) + epszL*trOC1(th,B(1,1),B(3,3))))./...
    ((trOC1(th,A(1,1),A(3,3))*epszL).^(5/2)).*dS1(th), ...
    0,acos(1-2*k));
LzNL2=-3/(4)*integral(@(rho) ...
    k*rho.*(chiz*trOC2(rho,A(1,1),A(3,3)) + epszL*trOC2(rho,B(1,1),B(3,3)))./...
    (trOC2(rho,A(1,1),A(3,3))*epszL).^(5/2),...
    0,2*sqrt(k*(1-k)));
LzNL=LzNL1+LzNL2;
L_NL=diag([LtNL,LtNL,LzNL]);

D_L=gammaL*L_L*A;
D_NL=gammaNL*L_L*A + gammaL*(L_NL*A + L_L*B);

end