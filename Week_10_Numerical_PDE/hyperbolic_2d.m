% Note: This code is 99% AI-generated
% 2D wave equation (explicit centered differences)
clear; close all;

% Parameters
Lx = 1; Ly = 1;          % domain size
Nx = 101; Ny = 101;      % grid points
dx = Lx/(Nx-1); dy = Ly/(Ny-1);
x = linspace(0,Lx,Nx); y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);

c = 1.0;                 % wave speed
% CFL condition (2D)
dt_max = 1/( c*sqrt(1/dx^2 + 1/dy^2) );
dt = 0.9 * dt_max;       % safety factor
fprintf('dt = %.4e, dt_max = %.4e\n', dt, dt_max);

Tend = 10.0;
Nt = ceil(Tend/dt);

% Initial condition: Gaussian pulse in centre
u0 = exp(-200*((X-0.5*Lx).^2 + (Y-0.5*Ly).^2));
v0 = zeros(size(u0));    % initial velocity

% Optional source function handle s(x,y,t)
sfun = @(X,Y,t) 0*X;     % zero source; replace with @(X,Y,t) ... if needed

% Preallocate
u_prev = u0;             % u^{n-1} (initially assume same; first step handled separately)
u     = u0;              % u^{n}
u_new = zeros(size(u));

% Precompute convenience
cx = (c*dt/dx)^2;
cy = (c*dt/dy)^2;

% first time step: u^1
% compute laplacian of u0 (interior)
lap = zeros(Ny,Nx);
lap(2:end-1,2:end-1) = (u0(2:end-1,3:end) - 2*u0(2:end-1,2:end-1) + u0(2:end-1,1:end-2))/dx^2 ...
                     + (u0(3:end,2:end-1) - 2*u0(2:end-1,2:end-1) + u0(1:end-2,2:end-1))/dy^2;

u = u0 + dt.*v0 + 0.5*(c*dt)^2 .* lap + 0.5*dt^2 .* sfun(X,Y,0);

% Apply Dirichlet BCs (example: zero displacement on boundary)
u(:,1) = 0; u(:,end) = 0; u(1,:) = 0; u(end,:) = 0;
u_prev = u0;

% Visualization setup
figure('Units','normalized','Position',[.2 .2 .4 .5]);
h = surf(X,Y,u,'EdgeColor','none'); %axis equal tight;
colorbar; caxis([-0.2 0.2]); title('t = 0'); view([45 45]);
axis([0 1 0 1 -1 1]); set(gcf,'Color','w');
drawnow;

% Time stepping
for n = 1:Nt
    t = n*dt;
    % compute interior Laplacian
    lap(2:end-1,2:end-1) = (u(2:end-1,3:end) - 2*u(2:end-1,2:end-1) + u(2:end-1,1:end-2))/dx^2 ...
                         + (u(3:end,2:end-1) - 2*u(2:end-1,2:end-1) + u(1:end-2,2:end-1))/dy^2;

    % update interior points explicitly
    u_new(2:end-1,2:end-1) = 2*u(2:end-1,2:end-1) - u_prev(2:end-1,2:end-1) ...
        + (c*dt)^2 .* lap(2:end-1,2:end-1) + dt^2 .* sfun(X(2:end-1,2:end-1), Y(2:end-1,2:end-1), t);

    % boundary conditions (Dirichlet zero)
    u_new(:,1) = 0; u_new(:,end) = 0; u_new(1,:) = 0; u_new(end,:) = 0;

    % swap time levels
    u_prev = u;
    u = u_new;

    % plot occasionally
    if mod(n,round(Nt/200))==0
        set(h,'ZData',u);
        title(sprintf('t = %.4f', t));
        drawnow;
    end
end
