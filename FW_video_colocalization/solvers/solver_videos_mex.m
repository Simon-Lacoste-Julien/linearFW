function y = solver_videos_mex(x, E, id, vi)
% x : N vector
% E : contains the edges (nEdges x 2)
% id: Nvid x 2 matrix containing the begining and end index of each image

y   = zeros(size(x));
% for each video [parallelizable]
for i = 1 : size(id,1)
  y( id(i,1) : id(i,2) ) = solver_video_mex( x( id(i,1) : id(i,2)) , E{i}, vi( id(i,1) : id(i,2),2 ));
end

