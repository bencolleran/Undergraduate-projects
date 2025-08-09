function t = LS_gaussian(x,x0,fwhm,derivative)
    c=fwhm/(2*sqrt(2*log(2)));
    if derivative==0
        t=(1/(sqrt(2*pi)*c))*exp(-((x-x0).^2)/(2*c^2))
    elseif derivative==1
        t=-(x-x0).*(1/((c^3)*sqrt(2*pi))).*exp(-((x-x0).^2)/(2*c^2))
    else
        display('Invalid value of derivative')
    end
end