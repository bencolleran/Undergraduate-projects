function t = LS_lorentzian(x,x0,fwhm,derivative)
    gamma=fwhm/(sqrt(3));
    if derivative==0
        t=(2/(pi*(sqrt(3))*gamma))*((1+(4/3)*((x-x0)/gamma).^2)).^(-1)
    elseif derivative==1
        t=(-16.*(x-x0)./(3*pi*sqrt(3)*gamma^3)).*(1+(4/3)*((x-x0)/gamma).^2).^(-2)
    end
end