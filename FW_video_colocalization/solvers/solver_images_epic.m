function y = solver_images_epic(x, ids)
% x : N vector
% ids: Nimg x 2 matrix containing the begining and end index of each image
persistent idx 

if ~exist('idx','var')
  idx = ([1:(numel(y) / 20)]-1) * 20;
end

y   = zeros(size(x));

% fast reshape
[~, b] = min(reshape(x, 20, []));
y(b + idx) = 1;

%{
% for each video [parallelizable]
for i = 1 : size(ids,1)
  y( ids(i,1) : ids(i,2) ) = solver_image( x( ids(i,1) : ids(i,2) ) );
end
%}

