function [cloudtype prec] = conv_strat_separation(prec, w, qc, T, pm)
   % prec(:), w(:,:), qc(:,:), T(:,:), pm(:)
   % check out:
   % Sui, Tao, and Simpson 1994 \cite{sui1994tropical}  
   % Convective–stratiform rainfall separation by cloud content
   % cloudtype: convective=1 stratiform=0
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
   nx = numel(prec);
   idx = [nx,1:nx,1]; % deal with periodicity
   k600 = findplev(pm,60000);
   cloudtype = nan(1,nx);
   for ix=1:nx
     if any(qc(ix,:) > 0.005e-3) % 0.005 g/kg
       ic=1+ix; % the central grid point starts at second index of idx
       if isnan(cloudtype(ix)) % to avoid going through points that are already defined
         if ( prec(ix) > 20 )  % [mm/hr]
            cloudtype(ix) = 1; 
         else
           [~,imelt] = min(abs(T(ix,:)-273.15)); % melting level
           if ( prec(ix) > 0.1 ) % drile
             % In raining area, qc > 0.5 g kg^-1 or max w above 600hPa exceed 5m s^-1
             if ( any(qc(ix,:) > 0.5e-3) | any(w(ix,k600:end)>5) )
              cloudtype(ix) = 1;
             else
              cloudtype(ix) = 0;
             end
           elseif ( any(qc(ix,:) > 0.025e-3) | any(w(ix,1:imelt)>5) )
             % In non-raining area stratiform area qc > 0.025 g kg^-1 or max w exceed 5m s^-1 below melting level
             cloudtype(ix) = 1; 
           else
             cloudtype(ix) = 0;
           end
         end
       end
     end
   end
   % this is the loosest condition, therefore do it at the end to cover all possible conv points
   for ix=1:nx
     if any(qc(ix,:) > 0.005e-3) % 0.005 g/kg
       ic=1+ix; % the central grid point starts at second index of idx
       if ( prec(ix) >= 2*mean(prec(idx([ic-1,ic+1]))) ) 
             % if rainrateexceeds the average of the two neighboring grid points, then all three points
             % are convective% Sui use four neighboring for grid point 1.5km, we should just use two neighboring for 4km 
         cloudtype(idx([ic-1,ic+1])) = 1; 
       end
     end
   end
% Bruan's criteria is less rigorous
            % if local point precip hasn't reached surface
            % w > 3 m/s or qliq > 0.5 g/kg = 0.5e-3 kg/kg
%           if ( max(w)>3 | max(qc)>0.5e-3 )
%           elseif ( prec > 0.1 ) %[mm/hr] => 0.1mm/hr=2.4mm/day
