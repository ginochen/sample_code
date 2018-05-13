function [ill,xsquall,xlonsquall,ylatsquall,lonsquall,latsquall]=findsquall(lon,lat,ustorm,ixshift)
%  find the points with coordinate (i,j) that are along the storm normal vector
%  ustorm: storm velocity vector
%  ilonlat: index of lon and lat
%  lon, lat: the lons and lats of the entire x-y domain
%e=[-4,1];
syms x;  
%lon=360-[130:-4:110];lat=12:4:24;
xlon=(lon-lon(1))*110e3;  % convert deg to km 
ylat=(lat-lat(1))*110e3;
%vc=[xlon(ceil(numel(xlon)/2)),ylat(ceil(numel(ylat)/2))]; % storm center coord
%ixshift=-5; % shift the xlon center index
if ~exist('ixshift')
  ixshift=0;
end
vc=[xlon(ceil(numel(xlon)/2)+ixshift),ylat(ceil(numel(ylat)/2))]; % storm center coord
%gridr = (xlon(2)+ylat(2))/4;%
gridr = xlon(2)/2;%
%gridr = sqrt((xlon(2)/2)^2+(ylat(2)/2)^2)/2; % diaganol length of the grid box
ii=1;
for i=1:numel(lon); % lon index 
  for j=1:numel(lat); % lat index 
    vec = [xlon(i),ylat(j)]-vc; % vector relative to storm center
    d = norm(cross([ustorm',0],[vec,0]))/norm(ustorm); % shortest distance to mcs prog dir plane
    if d<=gridr  
      ilonlat(ii,:) = [i,j]; % lon lat pairs satisfying the closest cells to the plane  
      dis(ii) = d;
      ii=ii+1;
    end
%{
% solve the quadratic polynomial (x^2,x,1) root and see if it includes the diagnal width 
    dv = ([xlon(i),ylat(j)]-vc-ustorm*x);
    poly=dv(1)^2+dv(2)^2-gridr^2;
    range = double(solve(poly,x,'maxdegree',2));
%    range = double(solve(((lon(i)-vc(1)-e(1)*x)^2 + (lat(j)-vc(2)-e(2)*x)^2),...
%            x ,'maxdegree',2))
    
%    if (numel(range)==1 & gridr==range) || (gridr >= range(1) & gridr <= range(2))
    if all(isreal(range))% & numel(range)==2 
%    if range(1)==range(2) 
      ilonlat(ii,:) = [i,j]; 
      ii=ii+1;
    end
%}
  end
end
%scatter(lon(ilonlat(:,1)),lat(ilonlat(:,2)))
a = [ustorm',0];
b = [1,0,0]; % x-axis
if a(1)<=0
 ang = atan2d(norm(cross(-a,b)),dot(-a,b));
else
 ang = atan2d(norm(cross(a,b)),dot(a,b));
end
if ang<=45 % do pick merd
  ilontmp = unique(ilonlat(:,1)); % meridionally average
  for i = 1:numel(ilontmp)
    idx1 = find(ilonlat(:,1)==ilontmp(i));
    [~,idx2] = min(dis(idx1));
    ilatsquall(i) = ilonlat(idx1(idx2),2);
  end
  ilonsquall=ilontmp';
else % ang>45 do pick zonal
  ilattmp = unique(ilonlat(:,2)); % meridionally average
  for i = 1:numel(ilattmp)
    idx1 = find(ilonlat(:,2)==ilattmp(i));
    [~,idx2] = min(dis(idx1));
    ilonsquall(i) = ilonlat(idx1(idx2),1);
  end
  ilatsquall=ilattmp';
end
ill = [ilonsquall;ilatsquall]';
lonsquall=lon(ilonsquall);
latsquall=lat(ilatsquall);
xlonsquall=xlon(ilonsquall)-xlon(ilonsquall(1));
ylatsquall=ylat(ilatsquall)-ylat(ilatsquall(1));
xsquall=sqrt(xlonsquall.^2+ylatsquall.^2)/1000; % meters to km
xsquall=linspace(xsquall(1),xsquall(end),numel(xsquall));
scatter(lonsquall,latsquall);
xlim([lon(1),lon(end)])
ylim([lat(1),lat(end)]);



