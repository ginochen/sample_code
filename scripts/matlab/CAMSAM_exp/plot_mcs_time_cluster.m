%%%%%% plot after running trackMCS.m  %%%%%%%%%%%%%
  figure
  load mcs_clusters
%  [~, ~, ~, coastlat, coastlon] = plotmap(lat,lon,20,2,0);
%  scatter3(coastlon,coastlat,zeros(1,length(coastlat)),'k.'); 
  coast_centered(0); % center the coastline at 0 degrees
  hold on;
  cluster_cond=1;
  timefactor = 1/dt*3/24; % tranform to days
  axis equal % use actual aspect ratio
%%%%%%%%%%%%% Find Cluster with N-points %%%%%%%%%%%%%
  N1 = 3; N2 = 10
  for N = 3 % 2 steps is only 3hours apart, not 6 hours
    irg = find(nt4Cl>=N1 & nt4Cl<N2); % greater than or equal to N row indices
    irl = find(nt4Cl<N1); % less than N row indices
    if irg
      aa = mcslltC(nonzeros(mcslltClRowInd(irg,:)),1:3);
      scatter3(aa(:,1),aa(:,2),aa(:,3)*timefactor,'r.');
      aa = mcslltC(nonzeros(mcslltClRowInd(irl,:)),1:3);
      scatter3(aa(:,1),aa(:,2),aa(:,3)*timefactor,'k.');
    end
  end
  xlabel('longitude');ylabel('latitude');zlabel('time (days)')
  set(gcf,'color','w')
  ylim([-20 20])
  xlim([0 360])
  zlim([0 max(aa(:,3)*timefactor)])
  grid on; % grid minor <--- turns on minor grid
