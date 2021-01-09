function [M1n, p2p1n, rho2rho1n, T2T1n, p02p01n, p02p1n, M2n ] = NormShock(val, gamma, trigger)
%Spencer Moore
%normal shock Calculator 2021

%value is whatever you want to look up and type the trigger for what the
%value you want 1-6 being as shown in the outputs, not counting p02p1n

%test case
%{
clear;
clc;
val = 0.9582;
gamma = 1.4;
trigger = 3;
%}

syms M1 M2 p2p1 rho2rho1 T2T1 p02p01
assume(M1, {'real','positive'});
assume(M2, {'real','positive'});

eqnM2 = M2^2 == (1+((gamma-1)/2)*M1^2)/(gamma*M1^2 - (gamma-1)/2);
eqnrho = rho2rho1 == ((gamma+1)*M1^2)/(2+(gamma-1)*M1^2);
eqnp = p2p1 == 1 + (2*gamma)/(gamma+1)*(M1^2-1);
eqnT = T2T1 == (1 + (2*gamma)/(gamma+1)*(M1^2-1))*(((gamma+1)*M1^2)/(2+(gamma-1)*M1^2))^-1;
%theres probably a better way to do this but for the total pressure ratio
%we will use a modified entropy equation with the ideal gas constant
%divided out
eqnp0 = p02p01 == exp(-(gamma/(gamma-1))*log(T2T1)+log(p2p1));

switch trigger
    
    case 1

        M1n = val;
        
    case 2
        
        p2p1n = val;
        M1n = double(solve(subs(eqnp,p2p1,p2p1n),M1));
        
    case 3
        
        rho2rho1n = val;
        M1n = double(solve(subs(eqnrho,rho2rho1,rho2rho1n),M1));
        
    case 4
        
        T2T1n = val;
        M1n = double(solve(subs(eqnT,T2T1,T2T1n),M1));
        
    case 5
        
        p02p01n = val;
        M1n = double(solve(subs(subs(subs(eqnp0,T2T1,solve(eqnT,T2T1))...
            ,p2p1,solve(eqnp,p2p1)),p02p01,p02p01n),M1));
        
    case 6
        
        M2n = val;
        M1n = double(solve(subs(eqnM2,M2,M2n),M1));
        
end

M2n = double(solve(subs(eqnM2,M1,M1n),M2));
p2p1n = double(solve(subs(eqnp,M1,M1n),p2p1));
T2T1n = double(solve(subs(eqnT,M1,M1n),T2T1));
rho2rho1n = double(solve(subs(eqnrho,M1,M1n),rho2rho1));

[~ , p01p1] = IsenFlow(M1n, gamma, 1);
[~ , p02p2] = IsenFlow(M2n, gamma, 1);

%this is recalculated partially for case 5... oh well
p02p01n = p02p2*(p01p1)^-1*p2p1n;
p02p1n = p02p01n*p01p1;

end
