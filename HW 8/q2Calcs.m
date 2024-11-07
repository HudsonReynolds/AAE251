sympref('default');

syms c_d0 AR e m g rho S V x

cL = m * g / (rho * S * V^2);

cD = c_d0 + cL^2 / (pi * AR * e);

eqn = x == cL^(1/2) / cD;

solve(eqn, x);

displayFormula("eqn")

