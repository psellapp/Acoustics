function [St] = calc_rectangular_screech_tone(h,b,M_jet,M_design,a_inf,u_j,uc_uj)
% Calculate the screech tone of rectangular jet based on Tam(1988)
% 
% h,b - height and breadth of rectangular nozzle exit, [meters]
% u_j - fully-expanded jet exit velocity, [m/s]
% u_c - convective velocity, [m/s]
% a_inf - ambient speed of sound, [m/s] - value calculated at fully-expanded condition may not be a good fit. Using standard sea-level
% value 341 m/s seems to fit better with the measured data
% uc_uj - U_c/U_j
% M_jet, M_design - fully-expanded Mach number, design Mach number
% St - Strouhal frequency, [fh/u_j]
% All quantities with subscript:
%                               j - fully-expanded condition
%                               d - design condition
% 
% -Prabu, 12/2/2020

gamma = 1.4; %specific heat ratio for air
Mj_gamma = 1+(((gamma-1)/2)*(M_jet^2));
Md_gamma = 1+(((gamma-1)/2)*(M_design^2));

Aj_Ad = sqrt((M_design^2/M_jet^2)*((Mj_gamma/Md_gamma)^((gamma+1)/(gamma-1)))); %area ratio, A_jet/A_design
h_j = h*(1+((Aj_Ad-1)*(b/(h+b)))); %height of fully-expanded jet
b_j = b*(1+((Aj_Ad-1)*(h/(h+b)))); %breadth of fully-expanded jet

St = (uc_uj/(2*(1+((uc_uj*u_j)/a_inf))*sqrt(M_jet^2-1)))*sqrt((h_j/b_j)^2+1)/((((((Mj_gamma/Md_gamma)^((gamma+1)/(2*(gamma-1))))*(M_design/M_jet))-1)*(b/(b+h)))+1); %Strouhal number of screech frequency

end