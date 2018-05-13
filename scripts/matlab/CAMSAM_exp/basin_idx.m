%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose : Produce the global index set for the four tropical ocean basins
%            from using the variable LANDFRAC in nc outputs
%            (basins defined as Bretherton and Back 2004)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

landfrac = loadvar('LANDFRAC',camCase,1,0);
basin{1}='IND';  % indian ocean
cutoff_lat{1} = [[-20, 0],   [0,   20] ]; % lon=[40 to 120], lat=[-20 to 0], etc
cutoff_lon{1} = [[40,  120], [40,  100]];
basin{2}='WPAC'; % western pacific
cutoff_lat{2} = [[-20, 0],   [0,   20] ];
cutoff_lon{2} = [[120, 180], [100, 180]];
basin{3}='EPAC'; % eastern pacific
cutoff_lat{3} = [[-20, 10],  [10,  20] ];
cutoff_lon{3} = [[180, 300], [180, 271]];
basin{4}='ATL';  % atlantic
cutoff_lat{4} = [[-20, 10],  [-20, 10],  [10,  20]];
cutoff_lon{4} = [[  0, 20],  [290, 360], [271, 350]];  
basin{5}='ALL';
%
% find the ocean indices for all four basins
idx = [1:numel(landfrac)];
for ib = 1:numel(basin) % four basins
   iiocn{ib}=[];iinocn{ib}=[];
   for i = 1:numel(cutoff_lat{ib})/2
      iilon = find(lon>=cutoff_lon{ib}((i-1)*2+1) & lon<=cutoff_lon{ib}((i-1)*2+2));
      iilat = find(lat>=cutoff_lat{ib}((i-1)*2+1) & lat<=cutoff_lat{ib}((i-1)*2+2));
      basinlandfrac = landfrac(iilon,iilat);
      basinlandfrac(find(landfrac(iilon,iilat)==0))=NaN; % set the ocn to nan indices
      lfrac = landfrac;
      lfrac(iilon,iilat)=basinlandfrac;
      iiocn{ib} = cat(1,iiocn{ib},find(isnan(lfrac))); 
      if (i==numel(cutoff_lat{ib})/2)
         iinocn{ib} = setdiff(idx,iiocn{ib}); % the inverse index of iiocn (including land index and other basin index), 
                       %set the values for these indices to NaN, so we can plot the ones in the basin 
      end
   end
end
% obtain the indicies for a time-series of ocn basins 
for ib=1:numel(basin) 
   iiocnt{ib}=[]; 
   for it=1:numel(itimes); 
      iiocnt{ib} = cat(1,iiocnt{ib},iiocn{ib}+(it-1)*ngp); 
   end; 
end     
% save all the tropical ocean indices in one vector iiocnall
for ib = 1:numel(basin); iiocnall=cat(1,iiocnall,iiocn{ib}); end
iiocn{numel(basin)+1}=iiocnall;
save('basin_index.mat','iiocn','iinocn','basin')
save('basin_time-series_index','iiocnt','basin')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% test how the basin-wise variable field looks like
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%dv.prect=5;  varname{idv.prect} = 'PRECT';    var{idv.prect} = loadvar(varname{idv.prect},camCase,spd*1000,1);% (mm/day) 
%precc=mean(var{idv.prect},3);
%precctmp = precc;
%precctmp(iinocn{1})=NaN; %define the non-basin indices to NaN value
%contourf(precctmp');
