ii=1
for i=200:300
 load(['mcs_cluster_xy_var_' num2str(i) '.mat'],'vo','Nt','parm')
 if max(max(vo.pva{end}{ceil(Nt(i)/2)}(:,:,parm.i3km)))>0.6e-6% halflife max pva
   ii=ii+1;
 else
   disp(['no PVa i = ' num2str(i)])
 end
end
