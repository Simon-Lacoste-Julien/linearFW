function [x_t,f_t,res] = FW(x_0, A, b, fun_optim, opts) 
% [x_t,f_t,res] = FW(x_0, S_0, alpha_0, A, b, fun_optim, opts) 
% -> res: objective tracking
%
% This code runs standard FW 
% (so no need to maintain active set).
% 
%
% objective is: 0.5 x' A x + b' x
%
% [x_t, f_t] = FW(x_0, S_0, alpha_0, A, fun_optim, opts) 
% x_0 : variable initialization 
% A : the matrix in the cost function of the QP
% fun_optim : the linear minimization oracle function to solve the related LP:
%       min_{z in M} <z, grad(x_t) >
% for images : take the min of the integer solutions (eg z=(1,0,0) )
% for videos : shortest path algorithm
%
% opts.Tmax  % max number of iterations
% opts.TOL   % tolerance for convergence (FW gap stopping criterion)
% opts.verbose % whether to print progress
%
% res.primal = fvalues;
% res.gap = gap_values;
% res.x_t = x_t;

it           = 1;
minf = 0;
minx = [];
% init:

x_t         = x_0;

% tracking results:

fvalues = [];
gap_values = [];

fprintf('running plain FW, for at most %d iterations\n', opts.Tmax);

cost_fun =@(y) .5 * y' *A * y + b' * y;
% optimization: 

while it <= opts.Tmax
  it = it + 1; 

  % gradient:
  grad        = A * x_t + b;

  % cost function:
  f_t = cost_fun(x_t);

  % towards direction search:
  s_FW     = fun_optim( grad ); % the linear minimization oracle
  d_FW     = s_FW - x_t;

  % duality gap:
  gap = - d_FW' * grad;

  fvalues(it-1) = f_t;
  gap_values(it-1) = gap;
  
  if opts.verbose
      fprintf('it = %d -  f = %f - gap=%f\n', it, f_t, gap);
  end

  if gap < opts.TOL
    fprintf('end of FW: reach small duality gap (gap=%f)\n', gap);
    break
  end 
  
  % construct direction

  d = d_FW; 
  max_step = 1;
  
  % line search:

  step = -  (grad' * d) / ( d' * A * d );
  step = max(0, min(step, max_step ));
  
  if step < eps
      % not a descent direction???
      fprintf('ERROR -- not descent direction???')
      keyboard
  end
  
  % FW step:
  x_t = x_t + step * d; 
 
end

res.primal = fvalues;
res.gap = gap_values;
res.x_t = x_t;

end
