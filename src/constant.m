% The file contains the constants
global l0;
global nm;
global h;
global hbar;
global m0;
global eV;
global C;
global e_cgs;

l0 = 10; % characteristic length scale, in the unit of nm

nm = 1e-9; % 1 nm = 10-9 m

h = 6.62607e-34; % plank constant, 6.62607015×10−34 J⋅Hz−1

hbar = 1.05457e-34; % h/(2*pi), 1.054571817...×10−34 J⋅s, https://en.wikipedia.org/wiki/Planck_constant

m0 = 9.10938e-31; % electron rest mass, 9.1093837015(28)×10−31 kg, https://en.wikipedia.org/wiki/Electron_mass

eV = 1.60217e-19; % 1eV = 1.60217662x10-19 J

e_cgs = 4.80320e-10; % electrical charge in the unit of cgs, https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units#Physical_constants_in_CGS_units

C = hbar^2 / (2*m0*(l0*nm)^2) / eV; % scaling factor for energies, in the unit of eV