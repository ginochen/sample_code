function [z dzm dzi] = p2z(rho,dpi,dpm,g,nz);
   % p to z for CESM outputs
   for iz = 1:nz
%      if (iz > 1)
         dzm(iz)  = dpm(iz)/(rho(iz)*g); % dpm = rho * g * dzm from surface to mid-point and so on to the next mid-point, 
                                          % add surface geopotential height phis/g, maybe unecessary since near zero over the ocean 
         dzi(iz) = dpi(iz)/(rho(iz)*g); % dpm(1) is surface pressure to 1st midpoint pressure
%      else
%         dzm(1) = dpm(1)/(rho(1)*g); 
%         dzi(1) = dpi(1)/(rho(1)*g); 
%      end
      z(iz) = sum(dzm(1:iz)); % height at each midpoint for CRM
   end;
