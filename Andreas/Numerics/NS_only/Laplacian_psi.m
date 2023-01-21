function [L_psi] = Laplacian_psi(nx,ny,dx,dy)
%LAPLACIAN_P Discretizes the Laplacian of p with finite volume method
%   

% Number of internal nodes
N = nx * ny;
vector_index = @(j,k) k + (j-1)*ny;

% Coefficient vectors for all cells
e = ones(N,1);  % auxiliary vector
A_E = e./dx^2;  
A_W = e./dx^2;
A_N = e./dy^2;
A_S = e./dy^2;
A_P = (-2/dx^2 -2/dy^2).*e ;


%% Boundary conditions

% Left boundary - homogeneous Dirichlet: p_W = -p_P
m = 1:ny;
A_P(m) = A_P(m) - A_W(m) ;
A_W(m) = NaN; % mark coefficients outside of the matrix

% Right boundary - homogeneous Dirichlet: p_E = -p_P
j = nx;
k = 1:ny;
m = vector_index(j,k);

A_P(m) = A_P(m) - A_E(m) ;
A_E(m) = NaN; % mark coefficients outside of the matrix

% Top boundary - homogeneous Dirichlet: p_N = -p_P
j = 1:nx;
k = ny;
m = vector_index(j,k);

A_P(m) = A_P(m) - A_N(m);
A_N(m) = 0; % Remove link to bottom boundary

% Bottom boundary - homogeneous Dirichlet: p_S = -p_P
j = 1:nx;
k = 1;
m = vector_index(j,k);

A_P(m) = A_P(m) - A_S(m);
A_S(m) = 0; % Remove link to top boundary


%% Assemble coefficients into a matrix

% Shift upper and lower diagonals due to the stupid syntax of spdiags 
% (see "doc spdiags")
A_W(1:N-ny) = A_W(ny+1:N);
A_S(1:N-1 ) = A_S( 2  :N);
A_N(2:N   ) = A_N(1:N-1 );
A_E(ny+1:N) = A_E(1:N-ny);

L_psi = spdiags([A_W,A_S,A_P,A_N,A_E], [-ny,-1,0,1,ny], N, N);


end

