%% Script to reproduce the Lasso experiment in the paper:
%
% On the Global Linear Convergence of Frank-Wolfe Optimization Variants, 
% Simon Lacoste-Julien and Martin Jaggi, NIPS 2015.
%
% It reproduces top of Figure 2 and takes a few seconds to run on my
% laptop.
%
% minimize ||A*x -b||^2 s.t. x is in l1-ball of radius r


clear 


%% random test case
d=200; % # rows in A
n=500; % dimension of w = # columns of A = # variables in x
t=25; % 2t = number of non-zeros in the true solution xstar
rng('default'); rng(42);      % you need to outcomment this line if running in octave
A=randn(d,n);
xstar=[ones(t,1); -ones(t,1); -zeros(n-2*t,1)];
noise=0.1*randn(d,1);
b=A*xstar+noise;
one_norm_of_xstar = norm(xstar,1)

r = 20; % the regularization constraint imposed on the l_1-norm

  

Aall = r*horzcat(A,-A); % all 2n atoms, over the simplex

%%
opts.Tmax  = 1000; % max number of iterations
opts.TOL   = 1e-8; % tolerance for convergence
opts.verbose = true;
opts.pureFW = 1; % to not use away steps...

%% running FW:
[x_t,f_t, resFW] = FW(Aall, b, opts);

%% running AFW:
opts.pureFW = 0;
[x_t,f_t, resAFW] = FW(Aall, b, opts);

%% running pairwise FW:
opts.pureFW = 0;
[x_t,f_t, resPairFW] = FW_pair(Aall, b, opts);

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


