function [D_L,D_NL]=DNL_UX_SpheroidHemi_INT(a,epsL,epsNL)
c=1;
% Dyadic L: nonlinear uniaxial medium: doubly truncated spheroid
% Components of epsL
epstL=epsL(1,1);
epszL=epsL(3,3);
% Components of epsNL
chit=epsNL(1,1);
chiz=epsNL(3,3);
% Linear and nonlinear components of inv(eps)=A+B|E|^2
A=diag(1./diag(epsL));
B=-A*epsNL*A;

if a<=1/sqrt(2)
    k=a*sqrt(1 - a^2);
else
    k=1/2;
end

% gamma
gammaL=epszL/epstL;
gammaNL=chiz/epstL-(epszL*chit/epstL^2);

trOC1=@(th,Ct,Cz) a^2*sin(th).^2*Ct + (cos(th)-k).^2*Cz; 
trOC2=@(rho,Ct,Cz) rho.^2*Ct + k^2*Cz;
% Normal vector:
%   > u = (sin(th)*cos(ph)/a,sin(th)*sin(ph)/a,cos(th)/c)
normU=@(th) sqrt(sin(th).^2/a^2 + cos(th).^2); % norm of normal vector
dS1=@(th) a*sin(th).*sqrt(sin(th).^2 + a^2*cos(th).^2); % surface element

LtL=1/(4) * integral(@(th) ... % Transverse component > NonLinear
    sin(th).^2.*dS1(th)./...
    (normU(th).*(trOC1(th,A(1,1),A(3,3))*epszL).^(3/2)), ...
    0,pi/2);
LzL1=1/(2) * integral(@(th) ...% Z component > NonLinear
    (cos(th)).*(cos(th)-k).*dS1(th)./...
    (normU(th).*(trOC1(th,A(1,1),A(3,3))*epszL).^(3/2)), ...
    0,pi/2);
LzL2=1/(2)*integral(@(rho) ...
    k*rho./...
    (trOC2(rho,A(1,1),A(3,3))*epszL).^(3/2),...
    0,a);
% (i)+2(ii) -- two plane surfaces (top and bottom)
LzL=LzL1+LzL2;

LtNL=-3/(8) * integral(@(th) ... % Transverse component > NonLinear
    sin(th).^2.*dS1(th).*(chiz*trOC1(th,A(1,1),A(3,3)) + epszL*trOC1(th,B(1,1),B(3,3)))./...
    (normU(th).*(trOC1(th,A(1,1),A(3,3))*epszL).^(5/2)), ...
    0,pi/2);
LzNL1=-3/(4) * integral(@(th) ...% Z component > NonLinear
    (cos(th).*(cos(th)-k).*dS1(th).*(chiz*trOC1(th,A(1,1),A(3,3)) + epszL*trOC1(th,B(1,1),B(3,3))))./...
    (normU(th).*(trOC1(th,A(1,1),A(3,3))*epszL).^(5/2)), ...
    0,pi/2);
LzNL2=-3/(4)*integral(@(rho) ...
    k*rho.*(chiz*trOC2(rho,A(1,1),A(3,3)) + epszL*trOC2(rho,B(1,1),B(3,3)))./...
    (trOC2(rho,A(1,1),A(3,3))*epszL).^(5/2),...
    0,a);
% (i)+2(ii) -- two plane surfaces (top and bottom)
LzNL=LzNL1+LzNL2;

L_L=diag([LtL,LtL,LzL]);
L_NL=diag([LtNL,LtNL,LzNL]);

% disp([trace(L_L)-1/gammaL,trace(L_NL)-(epszL*chit-epstL*chiz)/(gammaL*epstL)^2]);

D_L=gammaL*L_L*A;
D_NL=gammaNL*L_L*A + gammaL*(L_NL*A + L_L*B);
 end