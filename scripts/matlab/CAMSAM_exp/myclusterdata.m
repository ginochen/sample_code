function [C Z Y] = myclusterdata(X,Y,distance,method,criterion,cutoff,doecon)
% Purpose: Hierarchichal Clustering Algorithm, 
%  1. pairwise distance of row vectors of X by pdist()
%  2. horizontal dendrograms by linkage()
%  3. clustering by cluster()
% Caution: 
%  If OUT OF MEMORY, use dosingleEu to calc Euclidean distance with single
%  precision if precision is not that critical
% cutoff: cutoff clusters that are less than (not equal) the cutoff dist
%if ~dosavememory
  disp('start pair-wise distance calc')
  if doecon
    disp('use mydist()')
    Y = mydist(X,Y,distance); % takes 22min for 800e6 operations of cityblock
  else
    Y = pdist(X,distance);
  end
  disp('end pair-wise distance calc')
% distance :
% 'euclidean', 'squaredeuclidean', 'seuclidean',
% 'cityblock', 'minkowski', 'chebychev',
% 'mahalanobis', 'cosine', 'correlation',
% 'spearman', 'hamming', 'jaccard'
% method :
% 'average', 'centroid','complete','median',
% 'single', 'ward', 'weighted'
  disp('start linkage')
  Z = linkage(single(Y),method); % k-row = k-link := horiz lines in dendrogram
  disp('end linkage')
%{
else % not suitable for mcs tracking
  if ismember(method,{'ward','centroid','median'}) & strcmp(distance,'euclidean')
%   use 'ward', 'centroid','median' with 'euclidean' distance and set 'savememory' to 'on'
%   to avoid saving distance matrix Y, and use data matrix X directly
    disp('start linkage')
    Z = linkage(X,method,distance,'savememory','on');
    disp('end linkage')
  else
    error('use ''ward'',''centroid'',''median'' with ''euclidean''')
  end
end
%}

disp('start cluster')
C = cluster(Z,'criterion',criterion,'cutoff',cutoff);
disp('end cluster')
% criterion :
% 'inconsistent' or 'distance'
%{
function Y = mydist(X,Y,cutoff)
% single precision euclidean distance
tic
nx = size(X,1);
%cutoff = cutoff^2;
ii = 1;
%Y = zeros(1,nx*(nx-1)/2+1);
for i = 1:nx % loop over upper triangular matrix indices (i,j), calc all possible pairs of (lon,lat) distance
  for j = i+1:nx
    Y(ii) = sqrt(sum((X(i,:)-X(j,:)).^2));
    ii=ii+1;
  end
end
toc  
% XI: 1-by-n <--- single row of X
%}
