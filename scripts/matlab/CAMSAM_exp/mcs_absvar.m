function varo = mcs_absvar(vari)
% take abs value of vari and output as varo, where vari is a cell variable
for i=1:numel(vari); 
  if ~isempty(vari{i}); 
    varo{i} = cellfun(@abs,vari{i},'Un',0); 
  end; 
end
