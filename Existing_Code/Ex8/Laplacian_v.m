function [L_v,b_v] = Laplacian_v(Nx,Ny,BC,dx,dy)
%LAPLACIAN_V Discretizes the Laplacian of v with finite volume method
%   

% Number of internal nodes
ny = Ny-1;
N = Nx * ny;
vector_index = @(j,k) k + (j-1)*ny;

% Coefficient vectors for all cells
e = ones(N,1);  % auxiliary vector
A_E = e./dx^2;  
A_W = e./dx^2;
A_N = e./dy^2;
A_S = e./dy^2;
A_P = (-2/dx^2 -2/dy^2).*e ;



%% Boundary conditions

% Right-hand-side vector
b_v = zeros(N,1);


% Left boundary - requires interpolation (Ghost-point method)
m = 1:ny;

b_v(m) = b_v(m) - 2*A_W(m) .* BC.v.left;
A_P(m) = A_P(m) -   A_W(m) ; % v_W = 2*v_w - v_p

A_W(m) = NaN; % mark coefficients outside of the matrix


% Right boundary - requires interpolation (Ghost-point method)
j = Nx;
k = 1:ny;
m = vector_index(j,k);

b_v(m) = b_v(m) - 2*A_E(m) .* BC.v.right;
A_P(m) = A_P(m) -   A_E(m) ; % v_E = 2*v_e - v_P

A_E(m) = NaN; % mark coefficients outside of the matrix


% Top boundary 
j = 1:Nx;
k = ny;
m = vector_index(j,k);

b_v(m) = b_v(m) - A_N(m) .* BC.v.top(:);

A_N(m) = 0; % Remove link to bottom boundary


% Bottom boundary
j = 1:Nx;
k = 1;
m = vector_index(j,k);

b_v(m) = b_v(m) - A_S(m) .* BC.v.bottom(:);

A_S(m) = 0; % Remove link to top boundary



%% Assemble coefficients into a matrix

% Shift upper and lower diagonals due to the stupid syntax of spdiags 
% (see "doc spdiags")
A_W(1:N-ny) = A_W(ny+1:N);
A_S(1:N-1 ) = A_S( 2  :N);
A_N(2:N   ) = A_N(1:N-1 );
A_E(ny+1:N) = A_E(1:N-ny);

L_v = spdiags([A_W,A_S,A_P,A_N,A_E], [-ny,-1,0,1,ny], N, N);


end

