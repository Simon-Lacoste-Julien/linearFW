% Extra exploratory experiments for Appendix E (not in paper)
% just making sure that the theoretical bound is correct.

% suppose triangle with point at [-1,0],[0,0] and [cos(theta),sin(theta)]
% -> the smaller theta, the worse the width
% % width = sin(theta/2)
% diameter^2 = (1+cos(theta))^2+ sin(theta)^2 = 2*(1+cos(theta))
% or also: diameter = 2 cos(theta/2)
% so constant in rate is:
% width^2/diameter^2 = 1/4*(tan(theta/2))^2 
%
% Note: by looking instead at (<rhat, dhat>)^2 angle (tighter constant)
% we gain a factor of 4: it should just be (sin(theta/2))^2


% consider 1/2 ||x-x_0||^2 as objective with x_0 = [1/2,0];
% note: condition number of 1

xstar = [-0.5; 0];

grad = @(x) (x-xstar);
f = @(x) 0.5*norm(x-xstar)^2;

clear opts

opts.maxK = 500; 

opts.TOL = 1e-5;

%%
theta = pi/50; 

constant =  1/4*(tan(theta/2))^2;

V = [0,0; -1, 0; cos(theta), sin(theta)]'; % vertices for polytope

% drawing triangle:
figure
edges = [2,1; 1,3; 2,3];
for i = 1:3
    e = edges(i,:);
    plot(V(1,e), V(2,e),'-pb');
    hold on
end

%%

opts.pureFW = 1; rFW = AFW(V, grad, f, opts);
opts.pureFW = 0; rAFW = AFW(V, grad, f, opts);
rpFW = pFW(V, grad, f, opts);

%%
plot(rFW.X(1,1:50), rFW.X(2,1:50),'x--g')
if size(rAFW.X,2) < 50
    ind_end = size(rAFW.X,2);
else
    ind_end = 50;
end
plot(rAFW.X(1,1:ind_end), rAFW.X(2,1:ind_end),'x--r')

if size(rpFW.X,2) < 50
    ind_end = size(rpFW.X,2);
else
    ind_end = 50;
end
plot(rpFW.X(1,1:ind_end), rpFW.X(2,1:ind_end),'s:k')

%%
figure
semilogy(rFW.primal, 'g')
hold on
semilogy(rAFW.primal, 'r')
semilogy(rpFW.primal, 'k')
legend({'FW', 'AFW', 'pFW'})

%% comparing ratio of improvement vs. theory:

ratio = rpFW.primal(2:end) ./ rpFW.primal(1:(end-1));
ratio_theo = 1-constant;
figure
plot(ratio./ratio_theo)
title('for pFW: h_{k+1}/h_k divided by (1-\rho) [should be <= 1]')

% AFW:
ratio = rAFW.primal(2:end) ./ rAFW.primal(1:(end-1));
ratio_theo = 1-constant/4;
figure
plot(ratio./ratio_theo)
title('for AFW: h_{k+1}/h_{k} divided by (1-\rho/4) [should be <= 1]')
