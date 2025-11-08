f = @(x,y) [-2*y(1)-y(2)+exp(-x); y(1)];        % dy/dx = f(x,y)
h = 0.05; tol = 1e-5; x = 0:h:1;                % Step size; declare x
if x(end) < 1, x = [x 1]; end
 
alpha = 1; beta = 2;                            % Declare alpha, beta
theta = [1, 4];                                 % [theta_old, theta_new]
ft = @(t) RK4(f,x,[t;alpha],h,1)-beta;          % Solve roots of ft(t) = 0
yb1 = ft(theta(1));                             % Compute ft(theta_old)

while abs(theta(end)-theta(end-1)) > tol        % Run Secant method
    yb2 = ft(theta(end));
    theta(end+1) = theta(end) - (theta(end-1)...
        - theta(end))*yb2/(yb1 - yb2); %#ok
    yb1 = yb2;
end
 
y = RK4(f,x,[theta(end);alpha],h,0);            % Numerical solution, y(x)
g = @(x) 0.5*exp(-x).*(x.^2+(4*exp(1)-3)*x+2);  % True solution, y(x)
xtrue = linspace(x(1),x(end),200);              % Equally spaced x points
ytrue = g(xtrue);                               % Compute true solution
 
plot(x,y(2,:),'bo-'); hold on;                  % Plot numerical sol'n
plot(xtrue,ytrue,'k-','LineWidth',1.2);    	    % Plot true sol'n
legend(sprintf('Numerical Sol''n (h = %.2f,\\theta = %.2f)',...
    h,theta(end)),'True Solution','Location','northeast'); 
hold off; grid on;
 
function y = RK4(f,x,y0,h,last)
    y = repmat(y0,[1 length(x)]);           % Initialize y
    for j = 2:length(x)                     % Run RK 4 method
        k1 = h*f(x(j-1), y(:,j-1));
        k2 = h*f(x(j-1)+h/2, y(:,j-1)+k1/2);
        k3 = h*f(x(j-1)+h/2, y(:,j-1)+k2/2);
        k4 = h*f(x(j-1)+h, y(:,j-1)+k3);
        y(:,j) = y(:,j-1) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    if last == 1, y = y(2,end); end
end