% Declare the p, q, r functions
p = @(x) -2*ones(size(x));      % Declare p(x)
q = @(x) -1*ones(size(x));      % Declare q(x)
r = @(x) exp(-x);               % Declare r(x)
h = 0.01; x = 0:h:1;            % Step size; declare x
alpha = 1; beta = 2;            % Boundary conditions
if x(end) < 1, x = [x 1]; end

% Perform finite difference (AY = B)
A = diag(-(2+h^2*q(x))) + diag(1-p(x(1:end-1))*h/2,1) + ...
    diag(1+p(x(2:end))*h/2,-1);
B = (r(x)*h^2)'; B(1) = B(1) - (1+p(x(1))*h/2)*alpha;
B(end) = B(end) - (1-p(x(end))*h/2)*beta;
y = A\B;
 
g = @(x) 0.5*exp(-x).*(x.^2+(4*exp(1)-3)*x+2);  % True solution, y(x)
xtrue = linspace(x(1),x(end),200);              % Equally spaced x points
ytrue = g(xtrue);                               % Compute true solution
 
plot(x,y,'bo-'); hold on;                       % Plot numerical sol'n
plot(xtrue,ytrue,'k-','LineWidth',1.2);    % Plot true sol'n
legend(sprintf('Numerical Sol''n (h = %.3f)',h),...
    'True Solution','Location','southeast'); 
hold off; grid on; axis([0 1 0.5 2.5]);
