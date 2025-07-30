function [D_L,D_NL]=DNL_UX_SpheroidDT(a,k,epsL,epsNL)
% > a: equatorial radius
% > c: polar radius
% > k: truncation parameter
% > eps=epsL+epsNL*|E|^2
% > epsL: linear permittivity
% > epsNL: nonlinear permittivity
% > E2: |E|^2
% Dyadic L: nonlinear uniaxial medium: doubly truncated spheroid
% Components of epsL
% clear;
% epsL = diag([[1,1]*rand*5,rand*5]);
% epsNL=diag([[1,1]*rand*5,rand*5])*1e-4;
% a=rand*5;
% k=rand;

c=1;
epstL=epsL(1,1);
epszL=epsL(3,3);
% Components of epsNL
chit=epsNL(1,1);
chiz=epsNL(3,3);
% gamma
gammaL=epszL/epstL;
gammaNL=chiz/epstL-(epszL*chit/epstL^2);
% Linear and nonlinear components of inv(eps)=A+B|E|^2
A=diag(1./diag(epsL));
B=-A*epsNL*A;
% ======================
% Linear L: Using same expressions
Nu =1-a.^2.*gammaL;
Eta=k.^2-a.^2.*gammaL.*(k.^2-1);
    
LtL=(2*k.*Nu./(gammaL.*sqrt(Nu.*Eta))+...
    a.^2.*(log(-k.*sqrt(Nu)+sqrt(Eta))-log(k.*sqrt(Nu)+sqrt(Eta))))./(4*sqrt(Nu).^3);
    
LzL=1./gammaL-k./(gammaL.*sqrt(Eta)) +...
    (a.^2.*(-2*k.*sqrt(Nu.*Eta)-...
    Eta.*(log(-k.*sqrt(Nu)+sqrt(Eta))-log(k.*sqrt(Nu)+sqrt(Eta)))))./...
    (2*sqrt(Nu).^3.*Eta);
% ======================

LtNL = ((epszL * chit - epstL * chiz) / ...
(4*(gammaL*epstL)^2*(a^2*gammaL - 1)^(5/2)*(k^2 + a^2*gammaL - a^2*k^2*gammaL)^(3/2)))* ...
(sqrt(a^2*gammaL - 1)*k*...
(2*k^2 + 3*a^2*gammaL- 8*a^2*k^2*gammaL - 6*a^4*gammaL^2 + 6*a^4*k^2*gammaL^2) + ...
3*a^4*gammaL^2*(k^2 + a^2*gammaL - a^2*k^2*gammaL)^(3/2) * ...
atan(k*sqrt(a^2 * gammaL - 1)/sqrt(k^2 + a^2*gammaL - a^2*k^2*gammaL)));

LzNL=(epszL*chit-epstL*chiz)/(gammaL*epstL)^2 - 2*LtNL;
% ======================
L_L=diag([LtL,LtL,LzL]);
L_NL=diag([LtNL,LtNL,LzNL]);

% disp([trace(L_L)-1/gammaL,...
%     trace(L_NL)-(epszL*chit-epstL*chiz)/(gammaL*epstL)^2]);

D_L=gammaL*L_L*A;
D_NL=gammaNL*L_L*A + gammaL*(L_NL*A + L_L*B);
end