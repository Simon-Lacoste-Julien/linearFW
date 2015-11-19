%% Script to reproduce the flow QP experiment in the paper:
%
%   On the Global Linear Convergence of Frank-Wolfe Optimization Variants, 
%   Simon Lacoste-Julien and Martin Jaggi, NIPS 2015.
%
% It reproduces bottom of Figure 2 and takes less than 1 minute to run on my
% laptop.
%
% minimize 1/2 x' A x + b' x 
%      s.t. x belongs to some sort of flow polytope
%           arising in the video co-localization application
%           described in:
%  Efficient Image and Video Co-localization with Frank-Wolfe Algorithm,
%  Armand Joulin, Kevin Tang and Li Fei-Fei, ECCV 2014
%  -- their original code: http://ai.stanford.edu/~ajoulin/code/FW.zip


clear 

addpath solvers
addpath data

%%

is_imgs = false; % true : images / false : videos
is_test = true;

load('aeroplane_data_small.mat');

% get the LMO:
% linear function to optimize
% i.e. f(x) = argmin_{v in A} <v, x>
if is_imgs 
 [x_0, S_0, alpha_0, fun_optim, ids] = init_images(var_index);   
else
 [x_0, S_0, alpha_0, fun_optim, ids] = init_videos(var_index, edge_index,'mex');   
end


%%
opts.Tmax  = 2000; % max number of iteration
opts.TOL   = 1e-8; % tolerance for convergence
opts.verbose = true;

%% running FW:
[x_t,f_t, resFW] = FW(x_0, A, b,  fun_optim, opts);

%% running AFW:
[x_t,f_t, resAFW] = AFW(x_0, S_0, alpha_0, A, b,  fun_optim, opts);

%% running pairwise FW:
opts.pureFW = 0;
[x_t,f_t, resPairFW] = PFW(x_0, S_0, alpha_0, A, b,  fun_optim, opts);

%% plotting results:

figure
semilogy(resFW.gap,'b');
hold on
semilogy(resAFW.gap, 'k');
hold on
semilogy(resPairFW.gap, 'r');
legend({'FW', 'awayFW', 'pairFW'});


%%
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16)
xlabel('iteration','FontSize',20)
ylabel('gap','FontSize',20)

