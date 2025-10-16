f = @(x,y) [-x*y(1) - 3*y(2); y(1)];            % dy/dx = f(x,y)
h = 0.01; y0 = [0; 1]; x = 0:h:10;              % Step size, I.C., x
if x(end) < 10, x = [x 10]; end
 
y = repmat(y0,[1 length(x)]);           % Initialize y
for j = 2:length(x)                     % RK 4 method
    k1 = h*f(x(j-1), y(:,j-1));
    k2 = h*f(x(j-1)+h/2, y(:,j-1)+k1/2);
    k3 = h*f(x(j-1)+h/2, y(:,j-1)+k2/2);
    k4 = h*f(x(j-1)+h, y(:,j-1)+k3);
    y(:,j) = y(:,j-1) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
 
g = @(x) [exp(-x.^2/2).*(x.^3-3*x);             % True solution, y1(x)
          -exp(-x.^2/2).*(x.^2-1)];             % True solution, y2(x)
xtrue = linspace(x(1),x(end),200);              % Equally spaced x points
ytrue = g(xtrue);                               % Compute true solution

clf; subplot(211);
plot(x,y(1,:),'o-b'); hold on;                  % Plot numerical sol'n
plot(xtrue,ytrue(1,:),'k-','LineWidth',1.2);    % Plot true solution
hold off; axis([0 10 -2.5 1.5]);
subplot(212);
plot(x,y(2,:),'o-b'); hold on;                  % Plot numerical sol'n
plot(xtrue,ytrue(2,:),'k-','LineWidth',1.2);    % Plot true solution
hold off; axis([0 10 -1.5 2]);
legend(sprintf('Numerical Sol''n (h = %.2f)',h),...
    'True Solution','Location','northeast');
