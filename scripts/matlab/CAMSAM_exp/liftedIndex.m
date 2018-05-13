function [LI, tp500] = liftedIndex(Tp,pp,qv,p,t,tp500)
   % http://www.weathertap.com/guides/aviation/lifted-index-and-k-index-discussion.html
   % http://www.teachingboxes.org/avc/content/Severe_Weather_Indices.htm
   % https://www.ncl.ucar.edu/Support/talk_archives/2010/att-2526/sstats.f
   %  Tp : temp of the parcel
   %  pp : pres of the parcel [Pa]
   %  p  : vertical pressures on sigma-pressure levels [Pa] 
   %  t  : vertical temp on sigma-pressure levels
   %  qv : water vapor mixing ratio of the parcel
   %  lifting index LI = te500 - tp500
   %    Cp     = 1004; % [J kg^-1 K^-1]
   %    Rd     = 287.05; % [J kg^-1 K^-1]
   Rd   = 287;
   Cp   = 1004; %. or 1005.7 % specific heat dry air [J/kg/K]
   RCP  = Rd/Cp; 
   CPR  = Cp/Rd;
   P0   = 100000; % ref pressure for potential temp
   k500 = findplev(p,50000); % find the closest p-lev to 50000
   te500  = t(k500) + (t(k500+1)-t(k500))/(p(k500+1)-p(k500)) * (50000-p(k500)); % env temp interpolated to closest k500 level
%   if (nargin<6) % compare all parcels and find the most unstable (max(thelcl)) one and output it
   parm1  = (P0/pp)^RCP;
%      for ip = 1:length(qv) % number of parcels
   th     = Tp * parm1; % potential temp of the parcel 
   Td     = tdew(pp,qv); % dew point temp of the parcel
   tlcl   = t_lcl(Tp,Td); % temp at lcl 
   plcl   = P0 * (tlcl/th)^CPR; % pressure at lcl
   thelcl = theta_e(plcl,tlcl,qv,tlcl);
%      end
%      imaxp = find(thelcl == max(thelcl)); % max thelcl parcel index
%      thelcl = thelcl(imaxp(1)); % find the most unstable parcel, if more than one max, just use the first one since all values are the same
   tp500  = compT_fr_The(thelcl,50000); % temp of the parcel at p500
%   end
   LI  = te500 - tp500; % env - parcel 

function ip = findplev(p,plev); 
  % Purpose: find index of closest p to plev
  [~, ip] = min(abs(p-plev));

function Td = tdew(p,qv)
   % pressure in [Pa], temp in [K]
   RdRv = 0.622; 
   qv = qv+1e-8; 
   e = p*qv/(qv+RdRv);
   loge = log(e); 
   Td = (35.86*loge-4947.2325)/(loge-23.6837);


function Tlcl = t_lcl(Tp,Td)
   % The following code was based on Bolton (1980) eqn #15
   % and claims to have 0.1 K maximum error within -35 < T < 35 C
   %  Tp  = original parcel Temperature in Kelvin
   %  Td  = Temperature at Lifting Condensation Level (K)
   Tlcl = 1.0/(  1.0/(Td-56.0)  + log(Tp/Td)/800 )  + 56.0;


function th_e = theta_e(p, t, qv, tlcl)
   % from widepedia
   % https://en.wikipedia.org/wiki/Equivalent_potential_temperature
   %  e = p*qv/(qv+0.622);
   %  th_l = t*(100000/(p-e))^0.28*(t/tlcl)^(0.28*qv);
   %  th_e = th_l*exp((3036/tlcl-1.78)*qv*(1+0.448*qv));
   % from sstats.f
   qv    = qv + 1e-8;
   power = 0.2854*(1.0 - 0.28*qv);
   xx    = t * (100000.0/p)^power;
   p1    = 3.376/tlcl - 0.00254;
   p2    = (qv*1000.0) * (1.0 + 0.81*qv);
   th_e  = xx*exp(p1*p2);


function tp_p = compT_fr_The(thelcl,p)
   % compute parcel temp at pressure p using theta_e at lcl
   %  p: [Pa]
   %  thelcl: potential temp at LCL 
   %  tp_p: paracel temp at pressure p
   Tguess = (thelcl - 0.5 * max(thelcl-270, 0)^1.05)*(p/100000)^.2;
   epsilon=0.01;
   for iter=1:100
      w1 = r_sub_s(p,Tguess); % saturation mixing ratio
      w2 = r_sub_s(p,Tguess+1);
      tenu = theta_e(p,Tguess,w1,Tguess);
      tenup = theta_e(p,Tguess+1,w2,Tguess+1.);
      cor = (thelcl - tenu) / (tenup - tenu);
      Tguess = Tguess + cor;
      if ( (cor < epsilon) & (-cor < epsilon) ) 
         tp_p = Tguess; return
      end
   end
   thwlcl = theta_wetb(thelcl);
   tp_p = thwlcl*(p/100000.0)^0.286;

function qvs = r_sub_s(p,t)
   % this calls function e_sub_s which computes saturation
   % vapor pressure (Pa) and converts to sat. mixing ratio (kg/kg)
   %  p - pressure (pa)
   %  t  - temperature (k)
   %  qvs : staturation mixing ratio [kg/kg]
   RdRv = 0.622; 
   es = e_sub_s(t);
   qvs = RdRv*es/(p-es); 


function es = e_sub_s(t)
   % compute saturation vapor pressure (Pa) over liquid with
   % polynomial fit of goff-gratch (1946) formulation. (walko, 1991)
   c = [610.5851,44.40316,1.430341,.2641412e-1,.2995057e-3,.2031998e-5,.6936113e-8,.2564861e-11,-.3704404e-13];
   x = max(-80, t-273.16);
   es = c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)+x*(c(6)+x*(c(7)+x*(c(8)+x*c(9))))))));


function th_wetb = theta_wetb(thetae)
   % polynomial fit to data in  Smithsonian Meteorological Tables showing Theta-e and Theta-w
   c = [-1.00922292e-10, -1.47945344e-8, -1.7303757e-6, -0.00012709, 1.15849867e-6, -3.518296861e-9, 3.5741522e-12 ];
   d = [0.00000000,   -3.5223513e-10, -5.7250807e-8, -5.83975422e-6, 4.72445163e-8, -1.13402845e-10, 8.729580402e-14];
   x = min(475.0,thetae);
   if ( x <= 335.5 ) 
      th_wetb = 273.15 + c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)+x*(c(6)+ x*c(7) )))));
   else
      th_wetb = 273.15 + d(1)+x*(d(2)+x*(d(3)+x*(d(4)+x*(d(5)+x*(d(6)+ x*d(7) )))));
   end



