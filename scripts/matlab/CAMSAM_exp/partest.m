function partest
%parpool('local',2);
parpool('LSF1',32);
X = single([1,2,3;1,3,4;2,3,4;1,5,1;1,6,4;5,4,3;7,5,4;2,3,2;2,1,6;1,6,3]);
nx=size(X,1);
Y = zeros(1,nx*(nx-1)/2,'single');
  ii = 1;
  for i = 1:nx % loop over upper triangular matrix indices (i,j)
    for j = i+1:nx
      idx(ii,:) = [i,j];
      ii=ii+1;
    end
  end
 % for ii=1:size(idx,1)
  parfor ii=1:size(idx,1)
    Y(ii) = sum(abs(X(idx(ii,1),:)-X(idx(ii,2),:))); 
  end
delete(gcp('nocreate')); % close the parpool workers
