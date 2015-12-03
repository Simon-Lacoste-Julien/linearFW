function [x_t,f_t,res] = PFW(x_0, S_0, alpha_0, A, b, fun_optim, opts) 
% [x_t,f_t,res] = PFW(x_0, S_0, alpha_0, A, b, fun_optim, opts) 
% -> res: objective tracking
%
% This code runs pairwise FW (Algorithm 2 in paper). 
%
% objective is: 0.5 x' A x + b' x
%
% [x_t, f_t] = PFW(x_0, S_0, alpha_0, A, fun_optim, opts) 
% x_0 : variable initialization 
% S_0 : active set initialization 
% alpha_0 : weights initialization
%    x_0 is a dx1 vector (d = dimension)
%    alpha_0 is k x 1 vector of convex combination weights
%    we need:    x_0 = S_0 * alpha_0
%    S_0 is d x k, each column is an *atom* in the expansion for x_0
%        (it might also contains non-active atoms)
% Simplest initialization: x_0, alpha_0 = [1] and S_0 = x_0  -> this is
%       valid *only if* x_0 is an atom (because the hashing function
%                                       expects specific structure for
%                                       atoms)
%       By using a more general hashing function, you could also handle
%       an arbitrary point as an atom...
%
% A, b : defines the cost function for QP
% fun_optim : the linear minimization oracle function to solve the related LP:
%       min_{z in M} <z, grad(x_t) >
%   for images : take the min of the integer solutions (eg z=(1,0,0) )
%   for videos : shortest path algorithm
%
% opts.Tmax  % max number of iterations
% opts.TOL   % tolerance for convergence (FW gap stopping criterion)
% opts.verbose % whether to print progress
%
% res.primal = fvalues;
% res.gap = gap_values;
% res.number_away = number_away;
% res.number_drop = number_drop;
% res.S_t = S_t;
% res.alpha_t = alpha_t;
% res.x_t = x_t;
%
% IMPORTANT: to extend this code to a QP for another domain
% -- you just need to provide a different linear oracle (fun_optim)
% -- you also need to modify the hashing function below to the kind
%   of atoms appearing in your domain; it has to associate a unique
%   string (or number) to each possible atom which is output of your LMO.


it           = 1;
minf = 0;
minx = [];
% init:

x_t         = x_0;
S_t         = S_0;
alpha_t     = alpha_0;

% each column of S_t is a potential vertex
% I_active -> contains index of active vertices (same as alpha_t > 0)
mapping = containers.Map();
% this will map vertex hash to id in S_t (to see whether already seen vertex)
% alpha_t(i) == 0 implies vertex is not active anymore...
% alpha_t will be a weight vector so that x_t = S_t * alpha_t

% constructing mapping:
max_index = size(S_t,2); % keep track of size of S_t
for index = 1:max_index
    mapping(hashing(S_t(:,index))) = index; 
end
I_active = find(alpha_t > 0);

% tracking results:

fvalues = [];
gap_values = [];
number_drop = 0; % counting drop steps (max stepsize for away step)

fprintf('running pairwise FW, for at most %d iterations\n', opts.Tmax);

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
      fprintf('it = %d -  f = %g - gap=%g\n', it, f_t, gap);
  end

  if gap < opts.TOL
    fprintf('end of PFW: reach small duality gap (gap=%g)\n', gap);
    break
  end 
  
  % away direction search:

  if ~isempty(S_t) 
    id_A   = away_step(grad, S_t, I_active);
    v_A    = S_t(:, id_A);
    %d_A    = x_t - v_A;
    alpha_max = alpha_t(id_A);
  else
    fprintf('error: empty support set at step (it=%f)\n', it);
  end

  % construct pair direction (between towards and away):

    d = s_FW - v_A; %was for away: d_A;
    max_step = alpha_max; % was for away: alpha_max / (1 - alpha_max);
  
  % line search: (same as away, but different d)

  step = -  (grad' * d) / ( d' * A * d );
  step = max(0, min(step, max_step ));
  
  if abs(step - max_step) < 10*eps
      number_drop = number_drop+1;
  end

  if step < -eps
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
          %TODO: could possibly also remove it from S_t
      else
          alpha_t(id_A) = alpha_t(id_A) - step;
      end
      
   % towards part of the step:
      % alpha_t = (1-step)*alpha_t;
      
      % is this s_FW a new vertex?
      h = hashing(s_FW);
      if ~mapping.isKey(h)
          % we need to add vertex in S_t:
          max_index = max_index + 1;
          mapping(h) = max_index;
          S_t(:,max_index) = s_FW;
          id_FW = max_index;
          alpha_t(id_FW) = step; % this increase size of alpha_t btw
          I_active = [I_active, id_FW];
      else
          id_FW = mapping(h);
          if alpha_t(id_FW) < eps
              % we already had atom in 'correction poytope', but it was not
              % active, so now track it as active:
              I_active = [I_active, id_FW];
          end
          alpha_t(id_FW) = alpha_t(id_FW) + step;
      end
      
      % exceptional case: stepsize of 1, this collapses the active set!
      if step > 1-eps;
          I_active = [id_FW];
      end
  
  x_t = x_t + step * d; 
 
end


% returns the id of the active atom with the worst value w.r.t. the
% gradient
function id = away_step(grad, S, I_active)
    s = grad' * S(:,I_active);
    [~,id] = max(s);
    id = I_active(id(1));
end

% MODIFY THE FOLLOWING FUNCTION for the atoms in your domain
function h = hashing(sequence)
    % sequence is 0-1 vector
    % output is string with a for 0, b for 1.
    h = repmat('a',1,length(sequence));
    h(sequence == 1) = 'b';
end


res.primal = fvalues;
res.gap = gap_values;
res.number_drop = number_drop;
res.S_t = S_t;
res.alpha_t = alpha_t;
res.x_t = x_t;

end
