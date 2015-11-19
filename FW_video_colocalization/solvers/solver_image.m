function y = solver_image(x)

[~,id]   = min(x);
y        = zeros(size(x));
y(id(1)) = 1;
