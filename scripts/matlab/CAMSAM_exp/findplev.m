function ip = findplev(p,plev); 
  % Purpose: find index of closest p to plev
  [~, ip] = min(abs(p-plev));
