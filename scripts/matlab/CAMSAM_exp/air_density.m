
% Calculate Density of each gridbox
% Total density of air = rho = p/(R_d * T_rho) where T_rho ~~
% T(1+0.61qv - qc - qi)
%T = ncread([Case.dir Case.name],'CRM_T',start,count);
%qt = ncread([Case.dir Case.name],'CRM_QT',start,count); % total water (vapor + cloud liquid) mixing ratio
%qc = ncread([Case.dir Case.name],'CRM_QC',start,count);
%qi = ncread([Case.dir Case.name],'CRM_QI',start,count);
%qv = qt - qc - qi;
%T_rho = T.*(1+0.61*qv - qc - qi); % density temperature
%T_rho = reshape(T_rho,nx,nz);
%Rd = 287.058;
%for ip = 1:nz;
%   rho(:,ip) = crmlev(ip)./(Rd.*T_rho(:,ip));
%end
%% dz = dp / (rho * g) 
%dz(:,1) = (1000 - crmlev(1)) ./ (rho(:,1) .* g); %surface (~1000hPa) to first level thickness
%z(:,1) = dz(:,1);
%for ip = 2:nz
%   dz(:,ip) = (crmlev(ip-1) - crmlev(ip))  ./ (rho(:,ip) .* g);
%   z(:,ip)  = sum(dz(:,1:ip),2); % calculate the height of each grid point at constant pressure levels 
%end
