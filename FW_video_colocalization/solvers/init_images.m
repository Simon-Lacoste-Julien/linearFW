function [x_0, S_0, alpha_0, fun_optim, ids] = init_images(var_index)
    
N = size(var_index,1);
[~,ie] = unique(var_index(:,1:2), 'rows','last');
[~,ib] = unique(var_index(:,1:2), 'rows','first');
ids = [ib, ie];

fun_optim = @(x) solver_images(x, ids);

% initialization:
% WARNING the active set and their weights has to be intialize to
% such that x_0 =  S_0 * alpha_0';
n_imgs = size(ids,1);
boxes_per_img = N / n_imgs;
x_0     = ones(N,1) / boxes_per_img;
S_0     = zeros(N,boxes_per_img);
alpha_0 = ones(1, boxes_per_img) / boxes_per_img;
for i = 1 : boxes_per_img
    S_0(ids(:,1) + (i-1) * (N + 1)) = 1;
end
