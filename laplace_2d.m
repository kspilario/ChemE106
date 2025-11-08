% Original code by Sir Ric Roxas II
T_top = 100; T_bottom = 10; T_left = 75; T_right = 50;
N = 100;

xspan = linspace(0,1,N+1);
yspan = linspace(0,1,N+1);

m = length(xspan)-2;
n = length(yspan)-2; 
e = ones(m*n,1);

% Setup the A matrix
A = spdiags([e e -4*e e e],[-m -1 0 1 m],m*n,m*n);
ridx = (1:m)*n; ridx(end)=[]; A(ridx,ridx+1)=0;
lidx = (1:m)*n+1; lidx(end)=[]; A(lidx,lidx-1)=0;

% Setup the B vector
B = zeros(m*n,1);
B(1:m) = B(1:m)+T_top;
B(end-m+1:end) = B(end-m+1:end)+T_bottom;
B((1:m)*n) = B((1:m)*n)+T_right;
B(((1:m)-1)*n+1) = B(((1:m)-1)*n+1)+T_left;
B = -1*B;

% Solve the system
T = A\B;
T = reshape(T,m,n)';

% Plot results
imagesc(T); colormap('hot');
colorbar; axis off; box on;
set(gcf, 'Color', 'w');