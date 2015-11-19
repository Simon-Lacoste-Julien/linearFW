function res = AFW(V, grad, f, opts)
% function res = AFW(V, grad, f, opts);
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

naway = 0;
drop_finished = 0;

for k = 1:MAXK
  
  g_k = grad(x);
  
  % FW corner:
  [~, i_FW] = min(g_k'*V);
  s = V(:,i_FW);
  d_FW = s-x;
  
  % away corner:
  [~, i_A] = max(g_k'*V);
  v = V(:,i_A);
  d_A = x-v;
   
  gap = -g_k'*d_FW;
  g_A = -g_k'*d_A;
  
  X(:,k) = x;
  G(k) = gap;
  GA(k) = g_A;
  F(k) = f(x);
  
  if gap < TOL
      fprintf('CONVERGED in %d iterations.\n', k);
      break
  end
  
  if gap > g_A || o.pureFW
      %FW step:
      isFW = 1;
      d_k = d_FW;
      gamma_max = 1;
  else
      naway = naway+ 1;
      % away step:
      isFW = 0;
      d_k = d_A;
      gamma_max = alpha(i_A) / (1-alpha(i_A));
  end
  
  gap_dir = -g_k'*d_k;
  
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
  
  gamma = gamma_opt;
  x = x + gamma*d_k;
  
  if isFW
    alpha = (1-gamma)*alpha;
    alpha(i_FW) = alpha(i_FW) + gamma;
  else
    alpha = (1+gamma)*alpha;
    alpha(i_A) = alpha(i_A) - gamma;
  end
end

res.X = X;
res.gap = G;
res.primal = F;
res.gap_A = GA;
res.x = x;
res.alpha = alpha;
res.naway = naway;
res.drop_finished = drop_finished; % means early stopping...

end

function options = defaultOptions()

options = [];
options.pureFW = 0;
options.maxK = 500;
options.TOL = 1e-5;
options.seed = 0;

end % defaultOptions
