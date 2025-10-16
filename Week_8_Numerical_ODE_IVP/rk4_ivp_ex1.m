f = @(x) polyval([-2 12 -20 8.5],x);            % dy/dx = f(x,y)
h = 0.01; y0 = 1; x = 0:h:4;                    % Step size, I.C., x
if x(end) < 4, x = [x 4]; end
 
y = y0*ones(size(x));           % Initialize y
for j = 2:length(x)             % RK 4 method
    k1 = h*f(x(j-1));
    k2 = h*f(x(j-1)+h/2);
    k3 = h*f(x(j-1)+h/2);
    k4 = h*f(x(j-1)+h);
    y(j) = y(j-1) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
 
g = @(x) polyval([-0.5 4 -10 8.5 1],x);         % True solution, y(x)
xtrue = linspace(x(1),x(end),200);              % Equally spaced x points
ytrue = g(xtrue);                               % Compute true solution

clf; plot(x,y,'o-b'); hold on;                  % Plot numerical sol'n
plot(xtrue,ytrue,'k-','LineWidth',1.2);         % Plot true solution
hold off; axis([0 4 0 6]);
legend(sprintf('Numerical Sol''n (h = %.2f)',h),...
    'True Solution','Location','northwest');
