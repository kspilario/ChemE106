f = @(x,y) (x + 1)*sqrt(y);                     % dy/dx = f(x,y)
h = 0.01; y0 = 25; x = 5:h:10;                  % Step size, I.C., x
if x(end) < 10, x = [x 10]; end

y = y0*ones(size(x));            	            % Initialize y
for j = 2:length(x)
    y(j) = y(j-1) + f(x(j-1),y(j-1))*h;         % Euler's Method
end
 
g = @(x) ((x.^2 + 2*x - 15).^2)/16;             % True solution, y(x)
xtrue = linspace(x(1),x(end),200);              % Equally spaced x points
ytrue = g(xtrue);                               % Compute true solution
 
clf; plot(x,y,'o-b'); hold on;                  % Plot numerical sol'n
plot(xtrue,ytrue,'k-','LineWidth',1.2);         % Plot true solution
hold off; axis([5 10 0 700]);
legend(sprintf('Numerical Sol''n (h = %.2f)',h),...
    'True Solution','Location','northwest');
