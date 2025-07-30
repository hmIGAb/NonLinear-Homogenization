function [D_L,D_NL]=DNL_UX_Sphere(epsL,epsNL)
% clear; clc;
% epsL = diag([[1,1]*rand*5,rand*5]);
% epsNL=diag([[1,1]*rand*5,rand*5])*1e-4;

A=diag(1./diag(epsL));
B=-A*epsNL*A;

% Dyadic L: nonlinear uniaxial medium: doubly truncated spheroid
% Components of epsL
epstL=epsL(1,1);
epszL=epsL(3,3);
% Components of epsNL
chit=epsNL(1,1);
chiz=epsNL(3,3);

gammaL=epszL/epstL;
gammaNL=chiz/epstL-(epszL*chit/epstL^2);

if gammaL==1
    LtL=1/3;
    LzL=1/3;
    LtNL=2*(epszL*chit - epstL*chiz)/(5*epstL^2);
    LzNL=(epszL*chit - epstL*chiz)/(5*epstL^2);
else
    LtL = 1/(2*gammaL - 2*gammaL^2) + asec(sqrt(gammaL))/(2*(gammaL-1)^(3/2));
    LzL= 1/(gammaL - 1) - asec(sqrt(gammaL))/((gammaL-1)^(3/2));
    
    LtNL= (epstL*chiz - epszL*chit)*(sqrt(gammaL - 1)*(5*gammaL - 2) - 3*gammaL^2*asec(sqrt(gammaL)))/...
          (4*(gammaL - 1)^(5/2)*(gammaL*epstL)^2);
    LzNL= (epszL*chit - epstL*chiz)*(sqrt(gammaL - 1)*(2*gammaL + 1) - 3*gammaL  *asec(sqrt(gammaL)))/...
      (2*(gammaL - 1)^(5/2)*gammaL*epstL^2);
end
L_L = diag([LtL,LtL,LzL]);
L_NL = diag([LtNL,LtNL,LzNL]);

% disp([trace(L_L)-1/gammaL,trace(L_NL)-(epszL*chit-epstL*chiz)/(gammaL*epstL)^2]);

D_L=gammaL*L_L*A;
D_NL=gammaNL*L_L*A + gammaL*(L_NL*A + L_L*B);
end