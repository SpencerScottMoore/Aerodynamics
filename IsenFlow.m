function [Mn, p0pn, rho0rhon, T0Tn, AAstarn] = IsenFlow(val, gamma, trigger)
%Spencer Moore
%Isentropic Flow Calculator 2021

%test case

%{
clear;
clc;
val = 1.091;
gamma = 1.4;
trigger = 3;
%}

%value is whatever you want to look up and type the trigger for what the
%value you want 1-4 being as shown in the outputs, 5 being AAstar subsonic
%and 6 being supersonic

%symbolic and make assumptions consistant with the values
syms M p0p rho0rho T0T AAstar nu
assume(AAstar, 'positive');
assume(M, {'real','positive'});

%Define the isentropic equations
eqnp = p0p == (1+((gamma-1)/2)*M^2)^(gamma/(gamma-1));
eqnrho = rho0rho == (1+((gamma-1)/2)*M^2)^(1/(gamma-1));
eqnT = T0T == (1+((gamma-1)/2)*M^2);
eqnA = AAstar == sqrt((1/M^2)*((2/(gamma+1))*(1+((gamma-1)/2)*M^2))^((gamma+1)/(gamma-1)));
%eqnnu = nu

%known mach number
switch trigger
    case 1

    Mn = val;
    p0pn = double(solve(subs(eqnp,M,Mn),p0p));

    case 2
        
    p0pn = val;
    Mn = double(solve(subs(eqnp,p0p,p0pn),M));
    
    case 3
    
    rho0rhon = val;
    Mn = double(solve(subs(eqnrho,rho0rho,rho0rhon),M));
    
    case 4
    
    T0Tn = val;
    Mn = double(solve(subs(eqnT,T0T,T0Tn),M));
    
    case 5
    
    AAstarn = val;
    Mn = min(double(solve(subs(eqnA,AAstar,AAstarn),M)));
    
    case 6
    
    AAstarn = val;
    Mn = max(double(solve(subs(eqnA,AAstar,AAstarn),M)));
    
    otherwise 
    
    fprintf(['Incorrect trigger value, please use 1-6, for M, p0p, rho0rho,'...
        'T0T, AAstar subsonic, AAstar supersonic respectively']);
end

p0pn = double(solve(subs(eqnp,M,Mn),p0p));
rho0rhon = double(solve(subs(eqnrho,M,Mn),rho0rho));
T0Tn = double(solve(subs(eqnT,M,Mn),T0T));
AAstarn = double(solve(subs(eqnA,M,Mn),AAstar));
    
end

