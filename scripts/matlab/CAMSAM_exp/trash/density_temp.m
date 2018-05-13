function [rho T_rho RH] = density_temp(T,P,qT,qc,qi,compset)
% Purpose: use the ideal gas law to get the total density and the
% associated density tempature for a moist parcel at (T,P)

% P = 100000 [Pa] 
%Ra = 287; % [J K-1 kg-1]
%1/Ra = 0.0034843205
if (compset == 'm2005')
   qv = qT - qc;
   %   qvs = r_sub_s(P,T); % saturation mixing ratio, used for calculating RH = qv/qvs
   %end
elseif (compset == 'sam1mom')
   qv  = qT - qc - qi; % qt = total water; qc = cloud liquid; qi = cloud ice; qv = cloud vapor
end
if (nargout == 3) % even if T_rho is ~ not called it will be counted as nargout
   qvs = r_sub_s(P,T); % saturation mixing ratio, used for calculating RH = qv/qvs, r_sub_s is used in LI function
   RH = qv/qvs;
end
T_rho = T.*(1+0.61*qv - qc - qi); % density temperature
rho = P.*0.003484./T_rho ; % P = rho Ra T_rho => rho = P * 1/Ra * 1/T_rho

function qvs = r_sub_s(p,t)
%      this calls function e_sub_s which computes saturation
%      vapor pressure (Pa) and converts to sat. mixing ratio (kg/kg)
%          p - pressure (pa)
%          t  - temperature (k)
%          qvs : staturation mixing ratio [kg/kg]
   es = e_sub_s(t);
   qvs = 0.622*es/(p-es);

function es = e_sub_s(t)
%      compute saturation vapor pressure (Pa) over liquid with
%      polynomial fit of goff-gratch (1946) formulation. (walko, 1991)
   c = [610.5851,44.40316,1.430341,.2641412e-1,.2995057e-3,.2031998e-5,.6936113e-8,.2564861e-11,-.3704404e-13];
   x = max(-80, t-273.16);
   es = c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)+x*(c(6)+x*(c(7)+x*(c(8)+x*c(9))))))));

