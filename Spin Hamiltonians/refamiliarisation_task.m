A=g*bmagn/(2*boltzm);
x=linspace(0.1,300,300);
f = tanh(3*x*A);
%plot(x,f);

B=linspace(0,5,100);
P=tanh(A*B/100);
%plot(B,P)
gN=5.585694;
R=3E-10;
T=[-1,0,0;0,-1,0;0,0,2]*(mu0*gfree*gN*bmagn*nmagn)/(4*pi*R^3);
T(1,1)*1E-6/planck;
T(2,2)*1E-6/planck;
T(3,3)*1E-6/planck;

%largest value from T(3,3)
coupling=10E6*planck;
R=((mu0*gfree*gN*bmagn*nmagn)/(4*pi*coupling))^(1/3);
R*1E10;