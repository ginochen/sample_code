function [ilat_lim, lat_limits, h, coastlat,coastlon] = plotmap(lat,lon,latlim,maps,docbar,cbarPosition,icmap)
% plot the coast line in a limited latitude
  absdlat    = abs(abs(lat)-latlim); % |dlat|
  ilat_lim   = find(absdlat==min(absdlat));
  lat_limits = [lat(ilat_lim(1)) lat(ilat_lim(2))];
  switch maps 
  case 1 % mmap
     lon_limits = [lon(1) lon(end)];
     m_proj('Equidistant Cylindrical','lon',lon_limits,'lat',lat_limits);
     h = m_coast('color',[0 0 0],'linewidth',2.5);
     m_grid('box','on','fontsize',20);
     if (docbar)
     if nargin<5; cbarPosition='southoutside'; icmap=10; end
     if nargin<6; icmap=10; end
     colorbar(cbarPosition); 
     colormap(cmap(icmap)); % icmap=10
     end
  case 2 % matlab map
    lon_limits= [0 360]
  %%%%%%% matlab maps %%%%%%%%%%%%
  %axesm('miller','maplatlim',[-40 40],'maplonlim',[-20 60])
    projmeth='eqdcylin'; % 'eqdcylin': equidistant cylindrical; 'miller'
    axesm(projmeth,'maplatlim',lat_limits,'maplonlim',lon_limits)
    framem;
    gridm;
    mlabel;
    plabel
    load coastlines;
    h=plotm(coastlat,coastlon);
  end
