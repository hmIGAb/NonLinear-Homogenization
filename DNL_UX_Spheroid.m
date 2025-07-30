function[D_L,D_NL]=DNL_UX_Spheroid(a,c,epsL,chi)
% clear; clc;
% epsL = diag([[1,1]*rand*5,rand*5]);
% chi=diag([[1,1]*rand*5,rand*5])*1e-4;
% a=rand*5;
% c=rand*5;

epstL=epsL(1,1);
epszL=epsL(3,3);
chit=chi(1,1);
chiz=chi(3,3);
gammaL=epszL/epstL;
gammaNL=chiz/epstL-(epszL*chit/epstL^2);

A=diag(1./diag(epsL));
B=-A*chi*A;

LtL=c/2*(c/(gammaL*c^2-a^2*gammaL^2) + a^2*asec(a*sqrt(gammaL)/c)/(a^2*gammaL-c^2)^(3/2));
LzL=a^2*(sqrt(a^2*gammaL-c^2)-c*asec(a*sqrt(gammaL)/c))/(a^2*gammaL-c^2)^(3/2);

LtNL=-c*(-epszL*chit+epstL*chiz)*(c*(2*c^2-5*a^2*gammaL)*sqrt(a^2*gammaL-c^2)+3*a^4*gammaL^2*asec(a*sqrt(gammaL)/c))/...
    (4*gammaL^2*(a^2*gammaL-c^2)^(5/2)*epstL^2);
LzNL=-a^2*(-epszL*chit+epstL*chiz)*(sqrt(a^2*gammaL-c^2)*(c^2+2*a^2*gammaL)-3*a^2*c*gammaL*asec(a*sqrt(gammaL)/c))/...
    (2*gammaL*(a^2*gammaL-c^2)^(5/2)*epstL^2);

L_L=diag([LtL,LtL,LzL]);
L_NL=diag([LtNL,LtNL,LzNL]);

% disp([trace(L_L)-1/gammaL,...
%     trace(L_NL)-(epszL*chit-epstL*chiz)/(gammaL*epstL)^2]);

D_L=gammaL*L_L*A;
D_NL=gammaNL*L_L*A + gammaL*(L_NL*A + L_L*B);
end