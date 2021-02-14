clear all
syms x
warning off

s = 1650*1e-4; % contact seperation in cm
a = s/40; % radius of the contact in cm
d1 = 195*1e-4; % thickness of top layer (TL) in cm
d2 = 0.000001*1e-4; % thickness of interface layer
d3 = 0.1*1e-4; % thickness of bottom layer (Si) in cm

Rsh1 = 1e2; % Sheet resistance of the top layer, TL;
rho1 = Rsh1*d1; % Resistivity of the top layer in ohm.cm

% rhoc = 1e-3; % interface resistivity in ohm.cm2
rhoc = logspace(-4,4,100); % interface resistivity in ohm.cm2
rho2 = rhoc/d2; % Resistivity of the interface in ohm.cm
Rsh2 = rho2/d2; % Sheet resistance of the top layer;

Rsh3 = 0.12; % Sheet resistance of the bottom layer, Si
rho3 = Rsh3*d3; % Resistivity of the bottom layer, Si in ohm.cm

rho4 = 1e9; % Resistivity of the substrate, infinite

i = 1;
while (i <= length(rhoc))

    T4 = rho4; 
    w3 = @(x) (1 - exp(-2*d3*x))./(1 + exp(-2*d3*x));
    T3 = @(x) (w3(x)*rho3 + T4)./(1 + w3(x)*T4/rho3);
    w2 = @(x) (1 - exp(-2*d2*x))./(1 + exp(-2*d2*x));
    T2 = @(x) (w2(x)*rho2(i) + T3(x))./(1 + w2(x).*T3(x)/rho2(i));
    w1 = @(x) (1 - exp(-2*d1*x))./(1 + exp(-2*d1*x));
    T1 = @(x) (w1(x)*rho1 + T2(x))./(1 + w1(x).*T2(x)/rho1);
    A1 = @(x) T1(x)/rho1;

    R_int_s = integral(@(x) A1(x).*sin(x*a)./x.*(besselj(1,x*a)./(x*a) - besselj(0,x*s)/2),0,inf);
    R_int_2s = integral(@(x) A1(x).*sin(x*a)./x.*(besselj(1,x*a)./(x*a) - besselj(0,x*2*s)/2),0,inf);
    R_4PP(i) = pi/log(2)*(rho1/(2*a)*(4/pi)*(R_int_2s- R_int_s)); % Check Ehrenstein Pg. 32

    i = i + 1;
end

if length(rhoc)>1
    figure(1)
    semilogx(rhoc, R_4PP)
    hold all

    % dydx = 1./rhoc.*abs(gradient(rhoc)./gradient(R_4PP));
    dydx_2 = R_4PP(1:(length(rhoc)-1))./(rhoc(1:(length(rhoc)-1))).*abs(diff(rhoc)./diff(R_4PP));
    % dydx_3 = R_4PP.*abs(gradient(rhoc)./gradient(R_4PP));

    figure(3)
    loglog(rhoc(1:(length(rhoc)-1)),dydx_2)
    hold all
    % loglog(rhoc,dydx)

    % figure(3)
    % loglog(R_4PP,dydx_3)
    % hold all
    rhoc_2 = rhoc(1:(length(rhoc)-1));
    dydx_2 = dydx_2';
    rhoc = rhoc';
    R_4PP = R_4PP';
end
