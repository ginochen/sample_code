function Y = mydist(X,Y,distance)
% single precision euclidean distance
% TODO: try develop this code so don't need to save all points of X using cutoff distance
%      to save memory!!!! save the dist as logical/int8/sparse variable to save memory?
%      maybe do: Y = sparse(1,ny), and modify linkage() to deal with this!!
nx = size(X,1);
if strcmp(distance,'cityblock')
  disp('start cityblock')
  ii = 1;
%{
  for i = 1:nx % loop over upper triangular matrix indices (i,j)
    for j = i+1:nx
      idx(ii,:) = [i,j];
      ii=ii+1;
    end
  end
%  parfor ii=1:size(idx,1)
  for ii=1:size(idx,1)
    Y(ii) = sum(abs(X(idx(ii,1),:)-X(idx(ii,2),:))); 
  end
%}
  for i = 1:nx % loop over upper triangular matrix indices (i,j)
    for j = i+1:nx
      Y(ii) = sum(abs(X(i,:)-X(j,:)));
      ii=ii+1;
    end
  end
elseif strcmp(distance,'euclidean')
  disp('start euclidean')
  ii = 1;
  for i = 1:nx % loop over upper triangular matrix indices (i,j)
    for j = i+1:nx
      Y(ii) = sqrt(sum((X(i,:)-X(j,:)).^2));
      ii=ii+1;
    end
  end
end  
% XI: 1-by-n <--- single row of X
