function newipts = addpts(ilons,parm,N)
  % Add new lon around the existing MCS lon to fill the total N lon 
  % N: total # of longitudes to extract (only use even number)
  nlon=numel(parm.lon);
  Nh=ceil(N/2); % half of N
  newipts = zeros(1,N);
  N1=ceil(numel(ilons)/2); % # of 1st half of ilons
  N1c=Nh-N1; % remaining # 1st half of newilons
  N2=numel(ilons)-N1; % # of 2nd half of ilons
  N2c=Nh-N2; % remaining # 2st half of newilons
  newipts(Nh-N1+1:Nh+N2)=ilons; % add ilons to center of the N-element newilons vector
  for i=1:N1c % add to 1st half of newilons
    if ilons(1)-i>0
      newipts(N1c-i+1) = ilons(1)-i;
    else
      newipts(N1c-i+1) = nlon-i+1;
    end
  end
  for i=1:N2c % add to 2nd half of newilons
    if ilons(end)+i<=nlon
      newipts(Nh+N2+i) = ilons(end)+i;
    else
      newipts(Nh+N2+i) = ilons(end)+i-nlon;
    end
  end
