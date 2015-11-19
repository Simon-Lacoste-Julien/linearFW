function y = solver_videos(x, C, ids, var_index)
% x : N vector
% C is a Nvid cell containing sparse adjency matrix such tthat C(i,j) = 1 if the is an edge i->j
% ids: Nvid x 2 matrix containing the begining and end index of each image

y   = zeros(size(x));
% for each video [parallelizable]
for i = 1 : size(ids,1)
  y( ids(i,1) : ids(i,2) ) = solver_video( x( ids(i,1) : ids(i,2) ) , C{i}, var_index( ids(i,1) : ids(i,2),2 ));
end
