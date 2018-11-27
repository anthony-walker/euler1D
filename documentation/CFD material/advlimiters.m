%------------------------------------------------------------
% one-dimensional linear advection using flux limiters
% df/dt + U df/dx = 0; U = 1 and hence upwinded flux calculations are done
% Written By: Sourabh V. Apte (ME667)
% Several limiters are used and the equation advects the initial
% condition with sharp discontinuities using a "linear" advection 
% equation. The solution predicted by higher order schemes with limiters
% is compared to simple upwind scheme. 
%------------------------------------------------------------
n=101; nstep=300; length=1.0;h=length/(n-1);
dt=0.333333333*h; time=0.0; 
for i=1:n, 
    x(i)=h*(i-1);
end

f=zeros(n,1);fo=zeros(n,1);fh=zeros(n,1);foo=zeros(n,1);
fu=zeros(n,1);fuo=zeros(n,1);

%initial conditions...
for i=1:n, 
    if (x(i)>=0.0 & x(i)<=0.6), f(i)=exp(-100*(x(i)-0.3)^2); end; 
end;

for i=1:n, 
    if (x(i)>=0.6 & x(i)<=0.8), f(i)=1.0; end; 
end; 
fu=f;

hold on;plot(x,f,'k','linewidt',2);

% the big time loop...
for m=1:nstep,m,time
    hold off;plot(x,f,'linewidt',2); %axis([1 n 2.0, 4.5]); %plot solution
    hold on;plot(x,fu,'r','linewidt',2);axis([0 length 0 2]); pause(0.1)
    foo=f;
    for is=1:2 %two steps per timestep
        fo=f;
        for i=2:n-1, 
            bot=(fo(i)-fo(i-1));top=(fo(i+1)-fo(i)); r=top*bot/(bot^2+0.00001);
            % Limiters
             psi=max([0, min([2*r,0.5*(r+1),2])]); % monotonized central (MC)
            %psi=max([0, min([1.5*r,1]),min([r,1.5]) ]); % Sweby 
            % psi=max([0, min([2*r,1]),min([r,2])]); % superbee gtext('superbee')
            % psi=max([0, min([r,1])]); % minmod gtext('MINMOD')
            % psi=(r+abs(r))/(r+1); % van Leer
            % psi=0; % upwind
            % psi=r;
            %f_at_face = f_cell + psi*0.5*h*S; where S is slope at face
            % here kappa = -1 is used; i.e. only S- is used
            fh(i)=fo(i)+0.5*psi*(fo(i)-fo(i-1));
        end;
        bot=(fo(1)-fo(n-1));top=(fo(2)-fo(1)); r=top*bot/(bot^2+0.00001);
        psi=max([0, min([2*r,0.5*(r+1),2])]); % monotonized central (MC)
        %psi=max([0, min([1.5*r,1]),min([r,1.5]) ]); % Sweby 
        % psi=max([0, min([2*r,1]),min([r,2])]); % superbee
        % psi=max([0, min([r,1])]); % minmod
        % psi=(r+abs(r))/(r+1); %van Leer
        % psi=0; % upwind
        % psi=r;
        fh(1)=fo(1)+0.5*psi*(fo(1)-fo(n-1)); fh(n)=fh(1);
        for i=2:n-1,
            f(i)=fo(i)-(dt/h)*(fh(i)-fh(i-1) );
        end;
        f(1)=fo(1)-(dt/h)*(fh(1)-fh(n-1) );f(n)=f(1);
    end
    % average the solution for second oder time step
    for i=1:n,
        f(i)=0.5*(f(i)+foo(i)); 
    end
    
    %-----------------------
    % first order upwind
    fuo=fu; 
    for i=2:n-1,
        fu(i)=fuo(i)-(dt/h)*(fuo(i)-fuo(i-1)); 
    end
    fu(1)=fuo(1)-(dt/h)*(fuo(1)-fuo(n-1));fu(n)=fu(1);
    time=time+dt;
end;

% Replot initial conditions for comparison...
f=zeros(n,1);
for i=1:n, 
    if (x(i)>=0.0 & x(i)<=0.6), 
        f(i)=exp(-100*(x(i)-0.3)^2); 
    end; 
end; 
for i=1:n, 
    if (x(i)>=0.6 & x(i)<=0.8), 
        f(i)=1.0; 
    end; 
end; 
hold on;plot(x,f,'k','linewidt',2);

