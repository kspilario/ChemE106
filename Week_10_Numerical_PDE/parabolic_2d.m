% User-defined Parameters
alpha = 0.1;           % thermal diffusivity
tf = 3.0;              % Final time (sec)
L = 1.0;               % domain length (0 to L)
Nx = 20;               % no. of spatial endpoints           
Nt = 300;              % no. of time points

% Setup grid
x = linspace(0,L,Nx);  % grid for x
t = linspace(0,tf,Nt); % grid for t
dx = x(2)-x(1);        % spatial step size
dt = t(2)-t(1);        % time step size
r = alpha*dt/dx^2;     % lambda
close all;

% Stability check (FTCS requires r <= 0.5)
fprintf('r = %.3f (must be <= 0.5)\n', r);

% IC and BC (Dirichlet)
T = zeros(Nt, Nx);    % Initialize solution
T_initial = 0;
T(1,:) = T_initial;   % set initial condition
T(:,1) = 100;         % left boundary
T(:,end) = 50;        % right boundary

% FTCS update for one time step
for j = 2:Nt
    T(j, 2:end-1) = T(j-1, 2:end-1) + r*(T(j-1, 3:end)...
        - 2*T(j-1, 2:end-1) + T(j-1, 1:end-2));
end

% Plot results
figure; set(gcf, 'Color', 'w');
for j = 1:Nt
    plot(x, T(j,:), 'Color', 0.5*[1 1 1 1], ...
        'LineWidth',0.01); hold on;
    h = plot(x, T(j, :), 'k', 'LineWidth',1.5);
    title(sprintf('alpha = %.2f, Time = %.2f sec', ...
        alpha, t(j)));
    pause(0.1); delete(h);
end