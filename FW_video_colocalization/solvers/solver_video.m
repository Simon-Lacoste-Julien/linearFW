function y = solver_video(x, C, id_frame)
% x is a N vector representing each boxes in each frame of the video
% C is an NxN sparse adjency matrix such tthat C(i,j) = 1 if the is an edge i->j

N = size(x,1);
y = zeros(N,1);
% WARNING:
scores          = x;
argmin_score    = N;
parents         = zeros(N,1);

frame_end = id_frame(end);


% construct the paths
for i = 1 : N
  parents_of_i = find( C(:,i) );

  if ~isempty(parents_of_i) 
    [s, p] = min( scores(parents_of_i) );
    parents(i) = parents_of_i(p);
    scores(i) = scores(i) + s;

    % WARNING should take the argmin/argmin only on the last frame?
    if id_frame(i) == frame_end && (~exist('min_score','var') || scores(i) < min_score)
      min_score = scores(i);
      argmin_score = i;
    end
  end
end

% retrieve the best path:
while argmin_score ~= 0
  %fprintf('n = %d - score = %f\n', argmin_score, scores(argmin_score)); 
  y(argmin_score) = 1;
  argmin_score = parents(argmin_score);
end

