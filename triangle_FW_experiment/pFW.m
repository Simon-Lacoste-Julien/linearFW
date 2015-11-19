function res = pFW(V, grad, f, opts)
% function res = pFW(V, grad, f, opts);
% pairwise FW
%
% V(:,i) -> i^th vertex
%
% opts.pureFW=1 to get no away step (0 by default)
%  .maxK -> max number of iterations (500 by default)
%  .TOL -> stopping criterion on gap (1e-5 by default)
%  .seed -> seed for random starting point; set to 0 to use a deterministic
%  one; (0 by default)
%
% xstar = [-0.5; 0];
% 
% grad = @(x) (x-xstar);
% f = @(x) 0.5*norm(x-xstar)^2;

opt_default = defaultOptions();
if (nargin > 3)
    o = processOptions(opts, opt_default);
else
    o = opt_default;
end


%alpha0 = ones(size(V,2),1)/size(V,2); -> uniform starting point

if o.seed == 0
    alpha0 = ones(size(V,2),1);
    alpha0(end) = alpha0(end)*5;
    alpha0 = alpha0 / sum(alpha0);
else
    rng(o.seed)     % you need to change this line if running in octave
    alpha0 = rand(size(V,2),1);
    alpha0 = alpha0 / sum(alpha0);
end
x0 = V*alpha0;

MAXK = o.maxK; 

TOL = o.TOL;

x = x0;
alpha = alpha0;

drop_finished = 0;

for k = 1:MAXK
  
  g_k = grad(x);
  
  % FW corner:
  [~, i_FW] = min(g_k'*V);
  
  % away corner:
  [~, i_A] = max(g_k'*V);
  
  d_FW = V(:,i_FW)-x;
  d_k = V(:,i_FW) - V(:,i_A);
  
  gap = -g_k'*d_FW;
  gap_dir = -g_k'*d_k;
  
  X(:,k) = x;
  G(k) = gap;
  F(k) = f(x);
  
  gamma_max = alpha(i_A);
  
  gamma_opt = gap_dir / norm(d_k)^2;
  
  if gamma_opt >= gamma_max
      % drop step
%      gamma_opt = gamma_max;
%       fprintf('DROP STEP!\n');
%       k
%       keyboard
      % note that because we are on triangle, this drop step means we will
      % converge next step...
      fprintf('DROP step after %d iterations.\n', k);
      drop_finished = 1;
      break
  end
  
  if gap < TOL
      fprintf('CONVERGED in %d iterations.\n', k);
      break
  end
  
  gamma = gamma_opt;
  x = x + gamma*d_k;
  
  alpha(i_FW) = alpha(i_FW) + gamma;
  alpha(i_A) = alpha(i_A) - gamma;
  
end

res.X = X;
res.gap = G;
res.primal = F;
res.x = x;
res.alpha = alpha;
res.drop_finished = drop_finished; % means early stopping...

end

function options = defaultOptions()

options = [];
options.maxK = 500;
options.TOL = 1e-5;
options.seed = 0;

end % defaultOptions
