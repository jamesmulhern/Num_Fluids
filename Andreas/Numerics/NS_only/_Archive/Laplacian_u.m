function [L_u,b_u] = Laplacian_u(Nx,Ny,BC,dx,dy)
%LAPLACIAN_U Discretizes the Laplacian of u with finite volume method
%   


% Number of cells for u
nx = Nx - 1;
N  = nx * Ny;
vector_index = @(j,k) k + (j-1)*Ny;

% Coefficients
e = ones(N,1);
A_E = e./dx^2;  
A_W = e./dx^2;
A_N = e./dy^2;
A_S = e./dy^2;
A_P = (-2/dx^2 -2/dy^2).*e ;



%% Boundary conditions
% Right-hand-side vector
b_u = zeros(N,1);

% Left boundary
m = 1:Ny;
b_u(m) = b_u(m) - A_W(m) .* BC.u.left;
A_W(m) = NaN; % mark coefficients outside of the matrix

% Right boundary
j = nx;
k = 1:Ny;
m = vector_index(j,k);
b_u(m) = b_u(m) - A_E(m) .* BC.u.right;
A_E(m) = NaN; % mark coefficients outside of the matrix

% Top boundary - requires interpolation (Ghost-point method)
j = 1:nx;
k = Ny;
m = vector_index(j,k);

b_u(m) = b_u(m) - 2*A_N(m) .* BC.u.top(:);
A_P(m) = A_P(m) -   A_N(m) ;
A_N(m) = 0; % Remove link to bottom boundary

% Bottom boundary - requires interpolation (Ghost-point method)
j = 1:nx;
k = 1;
m = vector_index(j,k);

b_u(m) = b_u(m) - 2*A_S(m) .* BC.u.bottom(:);
A_P(m) = A_P(m) -   A_S(m) ;
A_S(m) = 0; % Remove link to top boundary



%% Assemble coefficients into a matrix

% Shift upper and lower diagonals due to the stupid syntax of spdiags 
% (see "doc spdiags")
A_W(1:N-Ny) = A_W(Ny+1:N);
A_S(1:N-1 ) = A_S( 2  :N);
A_N(2:N   ) = A_N(1:N-1 );
A_E(Ny+1:N) = A_E(1:N-Ny);

L_u = spdiags([A_W,A_S,A_P,A_N,A_E].*e, [-Ny,-1,0,1,Ny], N, N);


end

