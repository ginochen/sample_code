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
