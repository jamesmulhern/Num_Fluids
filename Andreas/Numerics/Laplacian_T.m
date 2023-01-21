function [L_T,b_T] = Laplacian_T(Nx,Ny,BC,dx,dy)
%Laplacian_T Discretizes the Laplacian of T with finite volume method
%   


% Number of cells for T
N = Nx * Ny;
vector_index = @(j,k) k + (j-1)*Ny;


% Coefficients
e = ones(N,1);
A_E = e./dx^2;  
A_W = e./dx^2;
A_N = e./dy^2;
A_S = e./dy^2;
A_P = (-2/dx^2 - 2/dy^2).*e ;



%% Boundary conditions
% Right-hand-side vector
b_T = zeros(N,1);

% Left boundary - requires interpolation (Ghost-point method)
m = 1:Ny;
b_T(m) = b_T(m) - 2 * A_W(m) .* BC.T.left;
% A_P(m)     = (3*1/dx^2+2*1/dy^2);
A_P(m) = A_P(m) -   A_W(m) ;
% A_W(m) = NaN; % mark coefficients outside of the matrix
A_W(m) = 0; % Remove link to rigth boundary


% Right boundary - requires interpolation (Ghost-point method)
j = Nx;
k = 1:Ny;
m = vector_index(j,k);
b_T(m) = b_T(m) - 2* A_E(m) .* BC.T.right;
% A_P(m)     = (3*1/dx^2+2*1/dy^2);
A_P(m) = A_P(m) -   A_E(m) ;
% A_E(m) = NaN; % mark coefficients outside of the matrix
A_E(m) = 0; % Remove link to left boundary


% Top boundary - requires interpolation (Ghost-point method)
j = 1:Nx;
k = Ny;
m = vector_index(j,k);

b_T(m) = b_T(m) - 2*A_N(m) .* BC.T.top(:);
% A_P(m)     = (2*1/dx^2+3*1/dy^2);
A_P(m) = A_P(m) -   A_N(m) ;
A_N(m) = 0; % Remove link to bottom boundary

% Bottom boundary - requires interpolation (Ghost-point method)
j = 1:Nx;
k = 1;
m = vector_index(j,k);

b_T(m) = b_T(m) - 2*A_S(m) .* BC.T.bottom(:);
% A_P(m) = 2*1/dx^2+3*1/dy^2;
A_P(m) = A_P(m) -   A_S(m) ;
A_S(m) = 0; % Remove link to top boundary

% Corners
%%
%south-west corner
j = 1;
k = 1;
m = vector_index(j,k);
% A_P(m)     = 3*(1/dx^2+1/dy^2);
A_P(m) = A_P(m) - A_S(m) - A_W(m);

% A(m,m)    = A_P(m);
% A(m,m+Ny) = AE(m);
% A(m,m+1)  = AN(m);
% b_T(m) = b_T(m) - (2*1/dx^2*BC.T.left(k) + 2*1/dy^2*BC.T.bottom(j));
b_T(m) = b_T(m) - (2*A_W(m)*BC.T.left(k) + 2*A_S(m)*BC.T.bottom(j));


%%
%south-east corner
j = Nx;
k = 1;
m = vector_index(j,k);
% A_P(m)     = 3*(1/dx^2+1/dy^2) ;
A_P(m) = A_P(m) - A_S(m) - A_E(m);

% A(m,m)    = A_P(m);
% A(m,m-Ny) = AW(m);
% A(m,m+1 ) = AN(m);
b_T(m) = b_T(m) - (2*A_E(m)*BC.T.right(k) + 2*A_S(m)*BC.T.bottom(j));

%%
%north-east corner
j = Nx;
k = Ny;
m = vector_index(j,k);
% A_P(m)     = 3*(1/dx^2+1/dy^2);
A_P(m) = A_P(m) - A_N(m) - A_E(m);
% A(m,m)    = A_P(m);
% A(m,m-Ny) = AW(m);
% A(m,m-1 ) = AS(m);
b_T(m) = b_T(m) - (2*A_E(m)*BC.T.right(k) + 2*A_N(m)*BC.T.top(j));

%%
%north-west corner
j = 1;
k = Ny;
m = vector_index(j,k);
% A_P(m)     = 3*(1/dx^2+1/dy^2);
A_P(m) = A_P(m) - A_N(m) - A_W(m);;
% A(m,m)    = A_P(m);
% A(m,m+Ny) = AE(m);
% A(m,m-1 ) = AS(m);
b_T(m) = b_T(m) - (2*A_N(m)*BC.T.left(k) + 2*A_W(m)*BC.T.top(j));



%% Assemble coefficients into a matrix

% Shift upper and lower diagonals due to the stupid syntax of spdiags 
% (see "doc spdiags")
A_W(1:N-Ny) = A_W(Ny+1:N);
A_S(1:N-1 ) = A_S( 2  :N);
A_N(2:N   ) = A_N(1:N-1 );
A_E(Ny+1:N) = A_E(1:N-Ny);

L_T = spdiags([A_W,A_S,A_P,A_N,A_E].*e, [-Ny,-1,0,1,Ny], N, N);


end

