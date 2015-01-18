clear all
dimensionlessPressure = 2;              %Pressure
nx=60;                        %Number of steps in space(x)
ny=60;                        %Number of steps in space(y)
T=.700;                        %Final time
dt=.001;                      %Width of time step
nt=ceil(T/dt);                %Number of time steps
dt=T/nt;                      %Corrected time step
L=2;                          %Length of pipe
H=1;                          %Height of pipe
cavityL=.3;                     %Cavity Height
cavityH=.6;                      %Cavity Length
dx=(L+cavityL)/(nx-1);                  %Width of space step(x)
dy=(H+cavityH)/(ny-1);                  %Width of space step(y)
x=0:dx:(L+cavityL);                     %Range of x(0,L) and specifying grid points
y=0:dy:(H+cavityH);                      %Range of y(0,H) and specifying grid points
u=zeros(ny,nx);             %Preallocating u
v=zeros(ny,nx);             %Preallocating v
p=zeros(ny,nx);               %Preallocating p
ustar=zeros(ny,nx);
vstar=zeros(ny,nx);
pstar=zeros(ny,nx);

Re=2000;
rho = 1;
nu = .1;

for it=0:nt

    bstar=zeros(ny,nx);
    for i=2:ny-1
        for j=2:nx-1
            if(i/ny*(H+cavityH)>cavityH || (j/nx*(L+cavityL) > ((L+cavityL)/2-cavityL/2) && j/nx*(L+cavityL) < ((L+cavityL)/2+cavityL/2)))
                partialVxx = (u(i, j+1)-u(i, j-1))/(2*dx);
                partialVyy = (v(i+1, j)-v(i-1, j))/(2*dy);
                partialVyx = (v(i, j+1)-v(i, j-1))/(2*dx);
                partialVxy = (u(i+1, j)-u(i-1, j))/(2*dy);
                bstar(i, j) = 1/dt*(partialVxx+partialVyy)-partialVxx*partialVxx+2*partialVxy*partialVyx+partialVyy*partialVyy;
            end
        end
    end

    for n=1:49
        pstar = p;
        for i=2:ny-1
            for j=2:nx-1
                if(i/ny*(H+cavityH)>cavityH || (j/nx*(L+cavityL) > ((L+cavityL)/2-cavityL/2) && j/nx*(L+cavityL) < ((L+cavityL)/2+cavityL/2)))
                    pstar(i, j) = (dx^2*dy^2*(rho*bstar(i, j))-dy^2*(p(i, j+1)+p(i, j-1))-dx^2*(p(i+1, j)+p(i-1, j)))/(-2*(dx^2+dy^2));
                end
            end
        end
        
        %Boundary Conditions:
        pstar(ceil(cavityH/(cavityH+H)*ny), 1:ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)) = pstar(ceil(cavityH/(cavityH+H)*ny)+1, 1:ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx));
        pstar(ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx):nx) = pstar(ceil(cavityH/(cavityH+H)*ny)+1, floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx):nx);
        pstar(1, ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx):floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)) = pstar(2, ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx):floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx));
        pstar(1:ceil(cavityH/(cavityH+H)*ny), ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)) = pstar(1:ceil(cavityH/(cavityH+H)*ny), ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)+1);
        pstar(1:ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)) = pstar(1:ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)-1);        
        pstar(:, 1) = 1;
        pstar(:, nx) = 0;
        pstar(ny, :) = pstar(ny-1, :);
        p = pstar;
    end

       
    for i=2:ny-1
        for j=2:nx-1
            if(i/ny*(H+cavityH)>cavityH || (j/nx*(L+cavityL) > ((L+cavityL)/2-cavityL/2) && j/nx*(L+cavityL) < ((L+cavityL)/2+cavityL/2)))
                secondDerivativeVxDx = (u(i, j+1)-2*u(i, j)+u(i, j-1))/(dx)^2;
                secondDerivativeVxDy = (u(i+1, j)-2*u(i, j)+u(i-1, j))/(dy)^2;
                Pdx = (p(i, j+1)-p(i, j-1))/(2*dx);
                VxVxDx = u(i, j)*(u(i, j+1)-u(i, j-1))/(2*dx);
                VyVxDy = v(i, j)*(u(i+1, j)-u(i-1, j))/(2*dy);
                ustar(i, j) = (-1/rho*Pdx+nu*(secondDerivativeVxDx+secondDerivativeVxDy)-VxVxDx-VyVxDy)*dt+u(i, j);
            
                secondDerivativeVyDx = (v(i, j+1)-2*v(i, j)+v(i, j-1))/(dx)^2;
                secondDerivativeVyDy = (v(i+1, j)-2*v(i, j)+v(i-1, j))/(dy)^2;
                Pdy = (p(i+1, j)-p(i-1, j))/(2*dy);
                VxVyDx = u(i, j)*(v(i, j+1)-v(i, j-1))/(2*dx);
                VyVyDy = v(i, j)*(v(i+1, j)-v(i-1, j))/(2*dy);
                vstar(i, j) = (-1/rho*Pdy+nu*(secondDerivativeVyDx+secondDerivativeVyDy)-VxVyDx-VyVyDy)*dt+v(i, j);
            end 
        end
    end
    u = ustar;
    v = vstar;
    u(ny, :) = 0; %Top Wall
    v(ny, :) = 0; %Top Wall
    u(1, ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx):floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)) = 0;
    v(1, ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx):floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)) = 0;
    u(ceil(cavityH/(cavityH+H)*ny), 1:ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)) = 0;
    v(ceil(cavityH/(cavityH+H)*ny), 1:ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)) = 0;
    u(ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx):nx) = 0;
    v(ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx):nx) = 0;
    
    u(1:ceil(cavityH/(cavityH+H)*ny), ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)) = 0; 
    v(1:ceil(cavityH/(cavityH+H)*ny), ceil(((L+cavityL)/2-cavityL/2)/(L+cavityL)*nx)) = 0; 
    u(1:ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)) = 0;
    v(1:ceil(cavityH/(cavityH+H)*ny), floor(((L+cavityL)/2+cavityL/2)/(L+cavityL)*nx)) = 0; 
end
uplot=u;
vplot= v;
Len = sqrt(uplot.^2+vplot.^2+eps);
uplot=uplot./Len; vplot=vplot./Len;
quiver(x,y,uplot,vplot,'k-')
%quiver(x,y,u,v,'k-')
axis equal
axis([0 (L+cavityL) 0 (H+cavityH)])
hold on 
pcolor(x,y,p)
colormap(jet)
colorbar
shading interp