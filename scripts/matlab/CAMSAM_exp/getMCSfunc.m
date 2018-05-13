function ip = findplev(p,plev); 
  % Purpose: find index of closest p to plev
  ap = abs(p-plev); % absolute difference of all p to plev
  ip = find(ap == min(ap)); % index of closest p  


function thetae = gettheta_e(p,T,qv)
  Td = tdew(p,qv);
  tlcl = t_lcl(T,Td);
  thetae = theta_e(p,T,qv,tlcl);

function [pi pm dpi dpm]  = hybrid2p(p0, ps,hyai,hybi,hyam,hybm,nz)
   % hybrid to p coordinate for CESM outputs
   for iz = 1:nz+1
      pi(iz) = hyai(iz)*p0 + hybi(iz)*ps; % interface pressure
   end
   for iz = 1:nz
      dpi(iz) = pi(iz) - pi(iz+1);
   end
   for iz = 1:nz
      pm(iz) = hyam(iz)*p0 + hybm(iz)*ps; % mid-point pressure
      if (iz > 1)
         dpm(iz) = pm(iz-1) - pm(iz);
      else
         dpm(1) = ps - pm(1);
      end
   end

function [z dzi] = p2z(rho,dpi,dpm,g,z0,nz);
   % p to z for CESM outputs
   for iz = 1:nz
      if (iz > 1)
         dzm(iz)  = dpm(iz)/(rho(iz)*g);
         dzi(iz) = dpi(iz)/(rho(iz)*g);
      else
         dzm(1) = dpm(1)/(rho(1)*g) + z0; % dpm = rho * g * dzm from surface to mid-point and so on to the next mid-point, 
                                          % add surface geopotential height phis/g, maybe unecessary since near zero over the ocean 
         dzi(1) = dpi(1)/(rho(1)*g) + z0;
      end
      z(iz) = sum(dzm(1:iz)); % height at each midpoint for CRM
   end;


function [mcs, C] = getMCSlonlat(qi_crm, qv_crm, prec, T, T_crm, lon, lat, ilatlim, p0, ps, hyai, hybi, hyam, hybm, nz, nx)
   % Get MCS cluster center lon-lat and ilon-ilat
   cond_plot=0;  % scatter plot for the MCS points and cluster center points
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Calc LI for all big iEOF-PC points %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       illct = 1; % counter for the lon lat indices with zonal lons over ocn only
       for ilon = 1:numel(lon)
         for ilat = ilatlim
           % decide if cloud-top temperature < 241K and if LS cloud fraction = 1
           [~, pm, ~, ~]  = hybrid2p(p0, ps(ilon,ilat), hyai, hybi, hyam, hybm, nz);
           %%%%%% find high cloud zonal points and the associated cloud-top height %%%%%%%%%%
           ictop=NaN(1,nx); tfcold=NaN(1,nx); % predefine no cold
           for ix = 1:nx
             tfcld = qi_crm(ilon,ilat,ix,:) > 1e-5; % true-false logic variable of cloud 1e-5kg/kg = 0.01g/kg 
             if ( any(tfcld) & T_crm(ilon,ilat,ix,max(find(tfcld)))<=260 ) % if high cloud
               ictop(ix) = max(find(tfcld)); % save cloud-top vertical index
               if (T_crm(ilon,ilat,ix,ictop(ix))<220) % high cloud with very cold top
                  tfcold(ix) = 1;
               end
             end
           end
           ixhcl  = find(~isnan(ictop)); % zonal index of high cloud   % Ex. suppose ixhcl = {1 4 5 7 8 9 31 32} => {{1},{4,5},{7,8,9},{31,32}} => {{1,31,32},{4,5},{7,8,9}}
           nxhcl = numel(ixhcl); % total zonal index of high cloud
           minPt = 10; % only if more than 12*4km=48km (approx pi*24^2 = 1.8e3km^2) of high cloud can be defined as HCS
           %%%%%% group the zonal points into high cloud systems (HCS) %%%%%%%%%%
           i = 1; j=1; % zonal index counter, HCS counter
           clear ixset
           if (nxhcl >= minPt) 
             while (i <= nxhcl)
               i1 = i;
               if (i~=nxhcl) % if the last index hasn't reached
                 while (ixhcl(i+1)-ixhcl(i)==1) % if neighbor
                   i=i+1; 
                   if (i==nxhcl); break; end % if the last index is reached
                 end
               end
               if ( i-i1+1 >= minPt ) 
                 ixset{j} = ixhcl(i1:i); % HCS set
                 j=j+1;
               end 
               i=i+1; 
             end
             if exist('ixset')
               nel = numel(ixset); % total # of HCSs
               %%%%%% merge the two boundary systems into one due to the periodic domain %%%%%%%%%%%%
               if ( any(ixset{1}==1) & any(ixset{end}==nx) ) % if 1 and 32 is in the set, then concat those two into one set
                 ixset{1} = [ixset{1},ixset{end}]; ixset{end}=[];
                 nel = nel - 1;
               end
               %%%%%% find at least one MCS out of the HCSs, stop looping once it's found %%%%%%%%%%
               for i = 1:nel 
                 ixhcs = ixset{i}; 
                 lrainpct = sum(prec(ilon,ilat,ixhcs)>1)/numel(ixhcs); %portion of HCS lightly raining
                 hrainpct = sum(prec(ilon,ilat,ixhcs)>6)/numel(ixhcs); %portion of HCS heavily raining
                 if ( lrainpct>=0.5 ) % more than 10*4km=40km of high cloud & more than 70% with rain > 1mm/hr
                   if ( any(tfcold(ixhcs)) ); % minimum cloud-top temp < 220
%                     if( any(prec(ilon,ilat,tfnnan)>6) ) % heavy rain > 6mm/hr exists
                     if ( hrainpct>=0.1 ) % more than 10% with heavy rain > 6mm/hr
                       mcs.lonlat(illct,1:2)  = [lon(ilon),lat(ilat)]; % remove all the overlapped lon-lat pairs
                       mcs.ilonlat(illct,1:2) = [ilon,ilat];
                       illct = illct+1; break % once counted, leave this LS grid point and go to the next
                     end
                   end    
                 end
               end
             end
           end
         end
       end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % test plot how the cluster looks like on a lon-lat map %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (cond_plot)
      figure;
      for i=1:size(C,1);
         scatter(mcs.lonlat(i,1),mcs.lonlat(i,2),'b'); hold on; 
         h(i)=text(mcs.lonlat(i,1),mcs.lonlat(i,2),num2str(C(i))); set(h(i),'FontSize',10) % label the clusters in number
      end; %pause
      scatter(mcs.llcentroids(:,1),mcs.llcentroids(:,2),'r'); hold on; pause
      % use the T cluster index as MCS cluster index
   end


function [var]       = load_reshape(Case, it, varname,  start, count, nlon, nlat, nx, nz); 
   var  = ncread([Case.dir Case.name{it}],varname,start,count); 
   var  = reshape(var,nlon,nlat,nx,nz);

function [div, vort] = div_vort(u, w, dx, dzi)
   [dudx dudz] = gradient(u,dx,dzi);
   [dwdx dwdz] = gradient(w,dx,dzi);
   div  = dudx + dwdz;
   vort = dudz - dwdx;


function [FFTke]     = fftke(u, w, nx, nzi)
   FFTu = fft2(u); % fft: columnwise zonal fft for each pressure level independently
   FFTw = fft2(w); % fft2: row and column-wise fft, i.e., fft(fft(uw))
   FFTmagu = abs([FFTu(1,:); 2*FFTu(2:nx/2,:); FFTu(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
   FFTmagw = abs([FFTw(1,:); 2*FFTw(2:nx/2,:); FFTw(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
   FFTmagu = ([FFTmagu(:,1), 2*FFTmagu(:,2:nzi/2), FFTmagu(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
   FFTmagw = ([FFTmagw(:,1), 2*FFTmagw(:,2:nzi/2), FFTmagw(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
   FFTke = 0.5*(FFTmagu.^2 + FFTmagw.^2);


function ilons     = getlongs(ilon, nlon, nlons_half, ilat, landfrac, noLand)
  % obtain longitude indices centered around ilon
   ilons = ilon-nlons_half:ilon+nlons_half;
   if (ilons(1)<=0);        ilons(find(ilons<=0))      = [nlon-sum(ilons<=0)+1:nlon]; end 
                 % if long fell out the left/right boundary, use the right/left boundary periodically
   if (ilons(end)>=nlon+1); ilons(find(ilons>=nlon+1)) = [1:sum(ilons>=nlon+1)];      end
   if (noLand) % if true, discard land points if landfrac is greater than zero
      if (sum(landfrac(ilons,ilat))~=0) % if there are land indices involved, discard this sample
         ilons = 0;
      end
   end


function [RH] = specific2relhum(T,P,q)
   qv = q/(1-q);
   qvs = r_sub_s(P,T);
   RH = qv/qvs*100;



function [rho T_rho RH ] = density_temp(T,P,qv,compset)
   % Purpose: use the ideal gas law to get the total density and the
   % associated density tempature for a moist parcel at (T,P)
   Rd   = 287;
   rRd  = 1/Rd;
   Rv   = 461.5; % [J K-1 kg-1]
   RdRv = Rd/Rv; % 0.622
   zvir = Rv/Rd - 1; %Rd/Rv = 461.5/287= 0.622, Rv/Rd=1.606 => zvir=0.61
   % P = 100000 [Pa] 
   %Rd = 287; % [J K-1 kg-1]
   %1/Rd = 0.0034843205
   if (nargout >= 3) % even if T_rho is ~ not called it will be counted as nargout
      qvs = r_sub_s(P,T); % saturation mixing ratio, used for calculating RH = qv/qvs, r_sub_s is used in LI function
      RH  = qv/qvs*100;
   end
   T_rho = T*(1+zvir*qv);
   %T_rho = T*(1+qv/0.622)/(1+qv+qc+qi);
   %T_rho1 = T.*(1+0.61*qv - qc - qi); % approximate density temperature
   rho = P.*rRd./T_rho ; % P = rho Ra T_rho => rho = P * 1/Ra * 1/T_rho (1/Rd=0.003484)


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
   if (nargin<6) % compare all parcels and find the most unstable (max(thelcl)) one and output it
      parm1  = (P0/pp)^RCP;
      for ip = 1:length(Tp) % number of parcels
         th     = Tp(ip) * parm1; % potential temp of the parcel 
         Td     = tdew(pp,qv(ip)); % dew point temp of the parcel
         tlcl   = t_lcl(Tp(ip),Td); % temp at lcl 
         plcl   = P0 * (tlcl/th)^CPR; % pressure at lcl
         thelcl(ip) = theta_e(plcl,tlcl,qv(ip),tlcl);
      end
      imaxp = find(thelcl == max(thelcl)); % max thelcl parcel index
      thelcl = thelcl(imaxp(1)); % find the most unstable parcel, if more than one max, just use the first one since all values are the same
      tp500  = compT_fr_The(thelcl,50000); % temp of the parcel at p500
   end
   LI  = te500 - tp500; % env - parcel 


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



function [cloudtype rainrate] = conv_strat_separation(prec, w, qc, qr, dbz, dbz_max, Ze, T, z, dzi)
   % cloudtype: conv=1 strat=0
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % convective-stratiform separation criteria by rainrate  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % convective criteria (tag the lat lon to see if overlapping events
   % occur) from Braun 2010 paper "Simulation and interpretation of the
   % genesis of Tropical Storm Gert (2005) as part of the NASA tropical
   % cloud systems and processes experiment":
   % if local point precip
   % > 20 mm/hr = 20/1000/3600 m/s = 5.6e-6 m/s <---- rain fall velocity
   % > 2*( averaged 24 neighboring precip rate ), then the 
   % point and its surrounding 24 points are convective (if its just 2D
   % then don't use this criteria 
   rainrate=prec;
   if ( prec > 20 ) % [mm/hr]
      cloudtype = 1; % convective
   else
      % if local point precip hasn't reached surface
      % w > 3 m/s or qliq > 0.5 g/kg = 0.5e-3 kg/kg
      if ( max(w)>3 | max(qc)>0.5e-3 )
         cloudtype = 1;
      elseif ( prec > 0.1 ) %[mm/hr] => 0.1mm/hr=2.4mm/day
         cloudtype = 0;
      else
         cloudtype = NaN;
      end
   end
      
   % stratiform:
   % if local point precip
   % > 0.1 mm/hr 
   % rainfall rate and reflectivity uses 
   % Z = a * R^b 
   % where Z [mm^6 m^-3] is the reflectivity, and R [mm h^-1] is the
   % rainfall rate, for Marshall-Palmer relation a=200, b=1.5
   % qc_thres = 0.0005; % [kg/kg]
   % w_thres = 3; %[m/s]
   % if ( qc > qc_thres) | w > w_thres )
   % RdRv = 0.622; 
   %   qv = qv+1e-8; 
   %   e = p*qv/(qv+RdRv);
   %   pd = p - e;
   %   rhoair = pd/(Rd*T);
   %   rhoh2o = qr*rhoair; %[kg/m3]

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % convective-stratiform separation criteria by radar reflectivity  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % http://www.wdtb.noaa.gov/courses/MRMS/ProductGuide/PrecipQualityControls/conv-strat-precip-sep.php
   % 1. max dBZ (<25km) >35 dBZ (since >45 is the higher threshold for convective, but
   %    actually only need 35)
   % 2. rainrate calc use Marshall-Palmer formula a=200, b=1.5
   %    https://en.wikipedia.org/wiki/DBZ_(meteorology)
   %    Z = a * R^b, dBZ \propto 10*log_10(Z/Z0)  
%   TC = T - 273.15; % temp in celcius
%   iz13 = find(abs(z-1300)  == min(abs(z-1300)));  % 1.3km height index
%   iz25 = find(abs(z-25000) == min(abs(z-25000)));
%   iz0C = find(abs(TC)      == min(abs(TC))); % 0C freezing level
%   % below 25km
%   if (any(qr(1:iz0C) > 0.01*1e-3)) % if rain mixing ratio > 0.01g/kg=0.01*1e-3kg/kg below freezing level
%      if (z(iz0C) < 2000) % freezing level < 2km
%         cloudtype=0; % stratiform
%      else % freezing level >= 2km
%         if (max(dbz(1:iz13)) > 45 | max(w(1:iz13))>3 | max(qc(1:iz13))>0.5*1e-3 ) 
%                                         % if max reflectivity below 1.3km > 45dBZ
%                                         % if max w below 1.3km > 3m/s 
%                                         % if max qc below 1.3km > 0.5g/kg=0.5*1e-3kg/kg 
%            cloudtype=1; % convective
%         else
%            cloudtype=0;
%         end
%      end
%   else % above 25km and below 250km (unfortunately spcam doesn't reach above 28km, so can't use this part)
%      if ( any(dbz(iz25:end)>0) )
%         vilwc = vil(Ze(iz25:end),dbz(iz25:end),dzi(iz25:end)); % vertically integrated liquid water content [kg/m^2]
%         vilwc = vil(Ze,dbz,dzi); % vertically integrated liquid water content [kg/m^2]
%      else
%         cloudtype=NaN;
%      end
%      if (exist('cloudtype')==0) % if cloudtype nonexist 
%         if (vilwc < 6.5)
%            cloudtype=0; % stratiform
%         else % vilwc >=6.5
%            izN10C = find(abs(TC+10) ==  min(abs(TC+10))) % -10C level
%            if ( dbz(izN10C) < 30 ) % if -10C level reflectivity < 30 dBZ
%               cloudtype=0; % stratiform
%            else % if -10C level reflectivity >= 30 dBZ 
%               cloudtype=1; % convective
%            end   
%         end
%      end
%   else
%      cloudtype=NaN;
%   end
%   if (dbz_max == -30)
%      rainrate = 0;
%   else
%      rainrate = (10^(dbz_max/10)/200)^(5/8); % [mm/hr] 8/5=1.6
%   end


   
function vilwc = vil(Ze,dbz,dzi)
   % Greene and Clark 1972
   % http://www.wdtb.noaa.gov/courses/MRMS/ProductGuide/SevereWeather/vertically-integrated-liquid.php
   % VIL:=vertically integrated liquid-water content [kg/m^2] derived
   % from radar reflectivity ( check consistency with actual LWC
   % Ze [mm^6 m^-3]
   % dbz = 10 * log10(Z/Z0)
   % VIL [kg m^-2]
%     dbz(dbz>56)=56; % to exclude contribuitions from ice, set to 56 dBZ
     ipos = find(dbz>0);
     Ze(dbz>56)=10^(56/10); % to exclude contribuitions from ice, set to 56 dBZ
     vilwc = 3.44*1e-6*dot(Ze(ipos).^(4/7),dzi(ipos));
%     vilwc = 3.44*1e-6*dot(dbz(ipos).^(4/7),dzi(ipos));

