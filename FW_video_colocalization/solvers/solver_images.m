function y = solver_images(x, ids)
% x : N vector
% ids: Nimg x 2 matrix containing the begining and end index of each image

y   = zeros(size(x));
% for each video [parallelizable]
for i = 1 : size(ids,1)
  y( ids(i,1) : ids(i,2) ) = solver_image( x( ids(i,1) : ids(i,2) ) );
end


