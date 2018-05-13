% Print figure;

% Print panel a;
% ----------------------------------------------------------------------------------->
% Define figure size;
fig_size = [19 9];

% Sets figure specifications;
fig = figure(1);
    fig.Color = [1 1 1];
    fig.Units = 'inches';
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 fig_size(1) fig_size(2)];
    fig.Position = [0 0 fig_size(1) fig_size(2)];
    fig.PaperPositionMode = 'auto';
    
% Define limits for the map;
lat_limits = [min(lat) max(lat)];
lon_limits = [min(lon) max(lon)];

% starts projection;
%m_proj('Miller Cylindrical','lon',lon_limits,'lat',lat_limits)
%m_proj('Mercator','lon',lon_limits,'lat',lat_limits)
m_proj('Equidistant Cylindrical','lon',lon_limits,'lat',lat_limits)

%m_pcolor(lon,lat,WW)
%caxis([-10 10])
%colormap(redblue3(100))

%shading flat
%m_coast('patch',[.9 .9 .9]);
m_coast('color',[0 0 0]);
m_grid('box','on','fontsize',4);

hold on
m_contourf(lon,lat,stdmap.etmp',30,'linestyle','none')
colorbar
colormap(cmap(0))

title('EABS standard deviation (t=1,700, p-ave)')
%cc = colorbar('southoutside','peer',gca);
%cpos = get(cc,'Position');
%   set(cc,'Position',[cpos(1) cpos(2)-.04 cpos(3) cpos(4)/2]);
%   set(get(cc,'ylabel'),'String', ['Energy transfer rate (\times 10^-^3 W/m^2) '],'fontsize',14);
%   set(cc,'fontsize',14)
