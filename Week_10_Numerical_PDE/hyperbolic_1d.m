% User-defined Parameters
c = 1.0;        % wave speed
tf = 2.0;       % Final time (sec)
L = 1.0;        % domain length (0 to L)
Nx = 100;       % No. of spatial endpoints
Nt = 200;       % No. of time points

x = linspace(0,L,Nx);   % grid for x
t = linspace(0,tf,Nt);  % grid for t
dx = x(2)-x(1);         % spatial step size
dt = t(2)-t(1);         % time step size
r = c*dt/dx;            % Courant number

% Stability check (CTCS requires r <= 1.0)
fprintf('r = %.3f (must be <= 1.0)\n', r);

% IC and BC (Dirichlet)
u = zeros(Nt+1, Nx);
u(1,:) = sin(4*pi*x);             % initial displacement (scenario 1)
%u(1,:) = cos(4*pi*x)-1;           % initial displacement (scenario 2)
%u(1,:) = exp(-((x-0.5)/0.05).^2); % initial displacement (scenario 3)
u(2,:) = u(1,:);                  % initial velocity = 0
u(:,1) = 0;                       % left boundary
u(:,end) = 0;                     % right boundary


% CTCS update for one time step
for j = 3:Nt+1
    u(j, 2:end-1) = 2*u(j-1, 2:end-1) - u(j-2, 2:end-1) + ...
        r^2*(u(j-1, 3:end) - 2*u(j-1, 2:end-1) + u(j-1, 1:end-2));
end

% Plot results
figure; set(gcf, 'Color', 'w');
for j = 1:Nt+1
    plot(x, u(j,:), 'Color', 0.5*[1 1 1 0.5], ...
        'LineWidth',0.01); hold on;
    h = plot(x, u(j, :), 'k', 'LineWidth',1.5);
    title(sprintf('c = %.2f, Time = %.2f sec', ...
        c, t(j)));
    axis([0 1 -2 2]);
    pause(0.1); delete(h);
end