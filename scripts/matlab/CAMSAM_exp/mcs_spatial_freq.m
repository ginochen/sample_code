function [mcs_spatialfreq, LON, LAT] = mcs_spatial_freq(mcsillt,t,lon,lat,startdate,enddate)
  % Purpose: Calc the total count of MCS at each spatial point in time (MAM, JJA, SON, DJF)
  % startdate: '0001-03-01-00000'
  % enddate:   '0001-05-31-75600'
  its = strmatch(startdate,t); % time index of start date in MAM
  ite = strmatch(enddate,t); % time index of end date in MAM
  irow = find(mcsillt(:,3)>=its & mcsillt(:,3)<=ite); % row indices that are in March April May
  ilons = min(mcsillt(:,1)):max(mcsillt(:,1));
  ilats = min(mcsillt(:,2)):max(mcsillt(:,2));
  ii=1;
  for i = 1:numel(ilons) 
    for j = 1:numel(ilats)
      mcs_spatialfreq(i,j) = sum(mcsillt(irow,1)==ilons(i) & mcsillt(irow,2)==ilats(j)); % [lon,lat,freq] for contour plot
      LON(i,j) = lon(ilons(i)); LAT(i,j) = lat(ilats(j));
    end
  end
   
