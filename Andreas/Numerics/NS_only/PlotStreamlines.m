function [psi] = PlotStreamlines(u,v,xcorners,ycorners)
%PLOTSTREAMLINES 
%   

dx = xcorners(2) - xcorners(1);
dy = ycorners(2) - ycorners(1);
nx = length(xcorners)-2;
ny = length(ycorners)-2;

rhs = diff(u,1,1)/dy - diff(v,1,2)/dx;
LHS = Laplacian_psi(nx,ny,dx,dy);

psi = LHS\rhs(:);

psi = [ zeros(1,nx+2) ;
        zeros(ny,1), reshape(psi,ny,nx), zeros(ny,1);
        zeros(1,nx+2) ];
    
figure();
contour(xcorners,ycorners,psi,'k');
xlabel('$x$',"Interpreter","latex");
ylabel('$y$',"Interpreter","latex");
title('Streamlines');

end

