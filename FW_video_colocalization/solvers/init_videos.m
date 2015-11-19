function [x_0, S_0, alpha_0, fun_optim, ids] = init_videos(var_index, C, type)

N = size(var_index,1);
[~,ie] = unique(var_index(:,1), 'rows','last');
[~,ib] = unique(var_index(:,1), 'rows','first');
ids = [ib, ie];

if ~exist('type','var') || ~strcmp(type, 'mex')
  fun_optim = @(x) solver_videos(x, C, ids, var_index);
else
  V = cell(size(C));
  for c = 1 : numel(C)
    [i,j,~] = find(C{c});
    V{c} = [i,j] - 1;
    if isempty(V{c})
        V{c} = -ones(1,2);
    end
    V{c} = int32(V{c});
  end
  ids = int32(ids);
  var_index = int32(var_index);
  fun_optim = @(x) solver_videos_mex(x, V, ids, var_index - 1);
end

%x_0 = fun_optim(-rand(N,1));
x_0 = fun_optim(-ones(N,1));

S_0 = x_0;
alpha_0 = 1;

