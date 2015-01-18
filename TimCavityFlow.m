clear all
dimensionlessPressure = 2;              %Pressure
nx=51;                        %Number of steps in space(x)
ny=51;                        %Number of steps in space(y)
T=.700;                        %Final time
dt=.001;                      %Width of time step
nt=ceil(T/dt);                %Number of time steps
dt=T/nt;                      %Corrected time step
L=1;                          %Length of cavity
H=1;                          %Height of cavity
dx=L/(nx-1);                  %Width of space step(x)
dy=H/(ny-1);                  %Width of space step(y)
x=0:dx:L;                     %Range of x(0,L) and specifying grid points
y=0:dy:H;                      %Range of y(0,H) and specifying grid points
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
            partialVxx = (u(i, j+1)-u(i, j-1))/(2*dx);
            partialVyy = (v(i+1, j)-v(i-1, j))/(2*dy);
            partialVyx = (v(i, j+1)-v(i, j-1))/(2*dx);
            partialVxy = (u(i+1, j)-u(i-1, j))/(2*dy);
            bstar(i, j) = 1/dt*(partialVxx+partialVyy)-partialVxx*partialVxx-2*partialVxy*partialVyx-partialVyy*partialVyy;
        end
    end

    for n=1:49
        pstar = p;
        for i=2:ny-1
            for j=2:nx-1
                pstar(i, j) = (dx^2*dy^2*(rho*bstar(i, j))-dy^2*(p(i, j+1)+p(i, j-1))-dx^2*(p(i+1, j)+p(i-1, j)))/(-2*(dx^2+dy^2));
            end
        end
        
        %Boundary Conditions:
        pstar(ny, :) = 0;
        pstar(1, :) = pstar(2, :); %dp/dy = 0 at y = 0
        pstar(:, 1) = pstar(:, 2); %dp/dx = 0 at x = 0
        pstar(:, nx) = pstar(:, nx-1); %dp/dx = 0 at x = 2 

        p = pstar;
    end

       
    for i=2:ny-1
        for j=2:nx-1
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
%            if(vstar(i, j) ~= 0)
%                i
%                j
%                Pdy
%            end
                
        end
    end
    u = ustar;
    v = vstar;
    u(:, 1) = 0; %Left wall
    u(:, nx) = 0; %Rightwall
    u(1, :) = 0; %Bottom wall
    u(ny, :) = 1; %Top Wall
    v(:, 1) = 0;
    v(:, nx) = 0;
    v(1, :) = 0;
    v(ny, :) = 0;
end
uplot=u;
vplot= v;
Len = sqrt(uplot.^2+vplot.^2+eps);
uplot=uplot./Len; vplot=vplot./Len;
quiver(x,y,uplot,vplot,'k-')
axis equal
axis([0 L 0 H])
hold on 
pcolor(x,y,p)
colormap(jet)
colorbar
shading interp

