function [x_t,f_t,res] = FW_pair(A, b, opts) 
% [x_t,f_t,res] = FW_pair(A, b, opts)
% -> res: objective tracking
%
% solves min ||x-b|| with x constrained to a finite atomic norm ball (here
% for the lasso, A is the set of scaled lasso columns and their negatives)
%
% This code runs pairwise FW (Algorithm 2 in paper).
%
% opts.Tmax  % max number of iterations
% opts.TOL   % tolerance for convergence (FW gap stopping criterion)
% opts.verbose % whether to print info

[d,n] = size(A); % n is the number of atoms here (twice the number of lasso data columns)

% init:
% alpha_t will be a weight vector so that x_t = S_t * alpha_t
alpha_t = zeros(n,1); alpha_t(1) = 1;
x_t         = A(:,1);  % this is tracking A\alpha in the lasso case


% each column of S_t = A(:,I_active) is a potential vertex
% I_active -> contains index of active vertices (same as alpha_t > 0)
% alpha_t(i) == 0 implies vertex is not active anymore...
% alpha_t will be a weight vector so that x_t = A * alpha_t
I_active = find(alpha_t > 0);

% tracking results:

fvalues = [];
gap_values = [];
number_drop = 0; % counting drop steps (max stepsize for away step)

fprintf('running pairwise FW, for at most %d iterations\n', opts.Tmax);

% optimization: 
it = 1;
while it <= opts.Tmax
  it = it + 1; 

  % cost function:
  f_t = objective_fun(x_t,b);
  % gradient = x-b:
  grad = grad_fun(x_t,b);

  % towards direction search:
  [id_FW, s_FW]   = LMO(A,grad); % the linear minimization oracle, returning an atom

  d_FW     = s_FW - x_t;

  % duality gap:
  gap = - d_FW' * grad;

  fvalues(it-1) = f_t;
  gap_values(it-1) = gap;

  if gap < opts.TOL
    fprintf('end of FW: reach small duality gap (gap=%f)\n', gap);
    break
  end 
  
  % away direction search:
  [id_A, v_A]   = away_step(grad, A, I_active);
  d_A    = x_t - v_A;
  alpha_max = alpha_t(id_A);

  % construct pair direction (between towards and away):

    d = s_FW - v_A; %was for away: d_A;
    max_step = alpha_max; % was for away: alpha_max / (1 - alpha_max);
  
  % line search: (same as away, but different d)
  %step = -  (grad' * d) / ( d' * A * d ); % was this for the video QP
  step = - (grad' * d) / ( d' * d );
  
    % simpler predefined stepsize
    %stepSize = 2 / (t+2);
  step = max(0, min(step, max_step ));
  
  if abs(step - max_step) < 10*eps
      number_drop = number_drop+1;
  end
  
  if opts.verbose
      fprintf('it = %d -  f = %f - gap=%f - stepsize=%f\n', it, f_t, gap, step);
  end
  
  if step < eps
      % not a descent direction???
      fprintf('ERROR -- not descent direction???')
      keyboard
  end
  
  % doing steps and updating active set:
  
  % away part of the step step:
      %alpha_t = (1+step)*alpha_t; % note that inactive should stay at 0;
      if abs(step - max_step) < 10*eps
          % drop step:
          number_drop = number_drop+1;
          alpha_t(id_A) = 0;
          I_active(I_active == id_A) = []; % remove from active set
      else
          alpha_t(id_A) = alpha_t(id_A) - step;
      end
      
   % towards part of the step:
      % alpha_t = (1-step)*alpha_t;
      
      alpha_t(id_FW) = alpha_t(id_FW) + step;
          
      I_active = find(alpha_t > 0); %TODO: could speed up by checking if id_FW is new
      
      % exceptional case: stepsize of 1, this collapses the active set!
      if step > 1-eps;
          I_active = [id_FW];
      end
  
  x_t = x_t + step * d; 
 
end

function obj = objective_fun(x,b)
    obj = sum((x-b).^2)/2;
    %obj = .5 * norm(x - b, 2); %gives a different result
end
function grad = grad_fun(x,b)
    grad = x - b;
end

%LMO returns the index of, and the atom of max inner product with the negative gradient
function [id, s] = LMO(A,grad)
    ips = grad' * A;
    [~,id] = min(ips);
    s = A(:,id);
end

% returns the id of the active atom with the worst value w.r.t. the
% gradient
    function [id, v] = away_step(grad, A, I_active)
    ips = grad' * A(:,I_active);
    [~,id] = max(ips);
    id = I_active(id(1));
    v = A(:,id);
end


res.primal = fvalues;
res.gap = gap_values;
res.number_drop = number_drop;
res.S_t = A(:,I_active);
res.alpha_t = alpha_t;
res.x_t = x_t;

end
