x=linspace(-5,5,1000);
%plot(x,LS_lorentzian(x,1,2/(pi),0),'r')
%hold on;
%plot(x,LS_gaussian(x,1,2*sqrt(2*log(2))/sqrt(2*pi),0),'b')
plot(x,LS_lorentzian(x,0,sqrt(3*sqrt(3)/(2*pi)),1),'r')
hold on;
plot(x,LS_gaussian(x,0,1.06,1),'b')