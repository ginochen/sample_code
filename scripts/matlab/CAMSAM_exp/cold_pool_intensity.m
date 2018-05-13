function [B buoy] = cold_pool_intensity(z,th_rho,th_rho_ave,parm)
  % calc cold pool buoyancy
  tmp = 9.8*(th_rho-th_rho_ave)./th_rho_ave; % [m s^-2] buoyancy by density pottemp, Feng et al 2015 Mechanism of convective cloud org...
  buoy = interp1(z, tmp, parm.zint,'linear');
  S = find(buoy<-0.003);  % height index set that are negative buoyancy
  if ~isempty(S) & S(1)==1 % the surface layer has negative buoyancy
    [s n] = continuousSet(S,length(S),0); % use only s{1} which is the continuous index set starting from the surface
    B = -2*sum(buoy(1:max(s{1})))*parm.dzi; % B [m/s]: cold pool strength (Rotunno et al 1988), 
                                              % max(s{1}): cold pool depth (the highest point that starts having negative buoyancy)
  else
    B=0; % not cold pool
  end
