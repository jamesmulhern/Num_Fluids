clearvars;
close all;
clc

% domain size
Lx = 1;
Ly = 1;
%
% number of cells
Nx = 10;
Ny = 10;
% number of lines
Nx_lines = Nx+1;
Ny_lines = Ny+1;
%
dx = Lx/Nx;
dy = Ly/Ny;
%
% cell corners
%
xcorners = 0:dx:Lx;
ycorners = 0:dy:Ly;
%
% cell centres
xcenter = 1/2*(xcorners(1:end-1)+xcorners(2:end));
ycenter = 1/2*(ycorners(1:end-1)+ycorners(2:end));
%
% velocity field at the corners

ux = xcorners * 0;
uy = -ycorners * 0;
%

% diffusivity coefficient
%
kappa = 0.01;
%
%Config
s = 1/10;
dt = s*dx^2/kappa;
t_end = 10;

% Boundary conditions
%
Tn = xcenter;
Ts = xcenter;
Te = ones(Nx, 1)';
Tw = xcenter * 0;

% Tn = xcenter;
% Ts = xcenter;
% Te = xcenter * 0;
% Tw = xcenter * 0;


%
% we are going to solve the problem in the form of a linear system
%
%      A T = B
%
% matrices initialization
%
M = Nx*Ny;
A = zeros(M,M);
B = zeros(M,1);
%

% these are help coefficients 
%
AP = zeros(M,1);
AE = AP;
AW = AP;
AS = AP;
AN = AP;

%
%%
%====================================
% Coefficients for interior cells
%====================================

for k = 1:Ny  
    for j = 1:Nx
        m = k + (j-1)*Ny;
        AW(m) = -0.5*dy*ux(j)/dx-kappa*dy/dx^2;
        AE(m) = 0.5*dy*ux(j+1)/dx -kappa*dy/dx^2;
        AS(m) = -0.5*dx*uy(k)/dy-kappa*dx/dy^2;
        AN(m) = 0.5*dx*uy(k+1)/dy -kappa*dx/dy^2;
        AP(m) = 2*kappa*(dx/dy^2+dy/dx^2);
    end
end
%%
% assign coefficients
% interior cells
for k = 2:Ny-1  
    for j = 2:Nx-1
        m = k + (j-1)*Ny;      
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1 ) = AS(m);
        A(m,m+1 ) = AN(m);
%         B(m)      = 0;
    end
end

%%
% fill in the coefficients for the boundary conditions
% south boundary
k = 1;
for j = 2:Nx-1
        m = k + (j-1)*Ny;  
        AP(m)     = kappa*(2*dy/dx^2+3*dx/dy^2);
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m+Ny) = AE(m);
        A(m,m+1)  = AN(m);
        B(m)      = 2*kappa*dx/dy^2*Ts(j);
end
%
%%
% west boundary
j = 1;
for k = 2:Ny-1
        m = k + (j-1)*Ny;
        AP(m)     = kappa*(3*dy/dx^2+2*dx/dy^2);
        A(m,m)    = AP(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1)  = AS(m);
        A(m,m+1)  = AN(m);
        B(m)      = 2*kappa*dy/dx^2*Tw(k);
end
%%
% north boundary
k = Ny;
for j = 2:Nx-1
        m = k + (j-1)*Ny;
        AP(m)     = kappa*(2*dy/dx^2+3*dx/dy^2);
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1)  = AS(m);
        B(m)      = 2*kappa*dx/dy^2*Tn(j);
end
%%
% east boundary
j = Nx;
for k = 2:Ny-1
        m = k + (j-1)*Ny;
        AP(m)     = kappa*(3*dy/dx^2+2*dx/dy^2);
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m-1)  = AS(m);
        A(m,m+1)  = AN(m);
        B(m)      = 2*kappa*dy/dx^2*Te(k);
end
%%
%south-west corner
j = 1;
k = 1;
m = k + (j-1)*Ny;
AP(m)     = 3*kappa*(dy/dx^2+dx/dy^2);
A(m,m)    = AP(m);
A(m,m+Ny) = AE(m);
A(m,m+1)  = AN(m);
B(m)      = 2*kappa*dy/dx^2*Tw(k) + 2*kappa*dx/dy^2*Ts(j);

%%
%south-east corner
j = Nx;
k = 1;
m = k + (j-1)*Ny;
AP(m)     = 3*kappa*(dy/dx^2+dx/dy^2) ;
A(m,m)    = AP(m);
A(m,m-Ny) = AW(m);
A(m,m+1 ) = AN(m);
B(m)      = 2*kappa*dy/dx^2*Te(k) + 2*kappa*dx/dy^2*Ts(j);

%%
%north-east corner
j = Nx;
k = Ny;
m = k + (j-1)*Ny;
AP(m)     = 3*kappa*(dy/dx^2+dx/dy^2);
A(m,m)    = AP(m);
A(m,m-Ny) = AW(m);
A(m,m-1 ) = AS(m);
B(m)      = 2*kappa*dy/dx^2*Te(k) + 2*kappa*dx/dy^2*Tn(j);

%%
%north-west corner
j = 1;
k = Ny;
m = k + (j-1)*Ny;
AP(m)     = 3*kappa*(dy/dx^2+dx/dy^2);
A(m,m)    = AP(m);
A(m,m+Ny) = AE(m);
A(m,m-1 ) = AS(m);
B(m)      = 2*kappa*dy/dx^2*Tw(k) + 2*kappa*dx/dy^2*Tn(j);

%%
% 
% now me make A a sparse matrix
A = sparse(A);
%

%% Create an L matrix for time dependant implicit solving
I = speye(M);
L = I - A * dt;

R = I + A * dt;
%% Loop

for t = dt:dt:t_end
    B = L\B;
%     Assign the BC again.
% fill in the coefficients for the boundary conditions
% south boundary
k = 1;
for j = 2:Nx-1
        m = k + (j-1)*Ny;  
        B(m)      = 2*kappa*dx/dy^2*Ts(j);
end
%

% west boundary
j = 1;
for k = 2:Ny-1
        m = k + (j-1)*Ny;
        AP(m)     = kappa*(3*dy/dx^2+2*dx/dy^2);
        A(m,m)    = AP(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1)  = AS(m);
        A(m,m+1)  = AN(m);
        B(m)      = 2*kappa*dy/dx^2*Tw(k);
end

% north boundary
k = Ny;
for j = 2:Nx-1
        m = k + (j-1)*Ny;
        B(m)      = 2*kappa*dx/dy^2*Tn(j);
end

% east boundary
j = Nx;
for k = 2:Ny-1
        m = k + (j-1)*Ny;
        B(m)      = 2*kappa*dy/dx^2*Te(k);
end

%south-west corner
j = 1;
k = 1;
m = k + (j-1)*Ny;
B(m)      = 2*kappa*dy/dx^2*Tw(k) + 2*kappa*dx/dy^2*Ts(j);


%south-east corner
j = Nx;
k = 1;
m = k + (j-1)*Ny;
B(m)      = 2*kappa*dy/dx^2*Te(k) + 2*kappa*dx/dy^2*Ts(j);


%north-east corner
j = Nx;
k = Ny;
m = k + (j-1)*Ny;
B(m)      = 2*kappa*dy/dx^2*Te(k) + 2*kappa*dx/dy^2*Tn(j);


%north-west corner
j = 1;
k = Ny;
m = k + (j-1)*Ny;
B(m)      = 2*kappa*dy/dx^2*Tw(k) + 2*kappa*dx/dy^2*Tn(j);
    t
end
% 
T = B;

%% we finally solve the linear system using backslash operator
% T = A\B;
%

%% we reshape
TM = reshape(T,Ny,Nx);

%% and plot
% [c, h] = contourm(xcenter,ycenter,TM);
% clegendm(c,h,4);

%% Surface plot
surf(xcenter,ycenter,TM);
colorbar;