S = 39.4

b = 19.78

AR = b^2 / S


powerMax = 1342.26 * 2 * 1000


[~,~,~,rho_h] = atmosisa(2709.062) 


[~,~,~,rho_0] = atmosisa(0) 


W = 10500

g = 9.81


e  = 0.7

m = 1

eta = 0.8


V = 156.464

cD0 = ((rho_h/rho_0)^m * P0 * eta - (2 * (W * g)^2) / (rho_h * S * V * pi * e * AR)) / (1/2 * rho_h * S * V^3)


