
%--------------------------------------------
% project: 1D Euler Equations: Lax Wendroff
% The Shock Tube Problem
% written by Sourabh V. Apte
% ME 667
%-------------------------------------------

nx=129; artvisc=0.5;
hold off
% gg=1.4;p_left=100000;p_right=1000;r_left=1;r_right=0.01;
gg=1.4;p_left=100000;p_right=10000;r_left=1;r_right=0.125;
xl=10.0;h=xl/(nx-1);time=0;
r=zeros(1,nx);ru=zeros(1,nx);rE=zeros(1,nx);p=zeros(1,nx);
rh=zeros(1,nx);ruh=zeros(1,nx);rEh=zeros(1,nx);ph=zeros(1,nx);

for i=1:nx,r(i)=r_right;ru(i)=0.0;rE(i)=p_right/(gg-1);end
for i=1:nx/2; r(i)=r_left; rE(i)=p_left/(gg-1); end

rh=r;ruh=ru;rEh=rE;ph=p;

dt=0.25*h/sqrt(1.4*max([700+p_right/r_right,700+p_left/r_left]) )

for istep=1:2000

for i=1:nx,p(i)=(gg-1)*(rE(i)-0.5*(ru(i)*ru (i)/r(i)));end

for i=2:nx-1					% prediction step
  rh(i)=0.5*(r(i)+r(i+1))-(0.5*dt/h)*(ru(i+1)-ru(i));
  ruh(i)=0.5*(ru(i)+ru(i+1))-...
	(0.5*dt/h)*((ru (i+1)^2/r(i+1))+p(i+1)-(ru (i)^2/r(i))-p(i));
  rEh(i)=0.5*(rE(i)+rE(i+1))-...
	(0.5*dt/h)*((rE(i+1)*ru (i+1)/r(i+1))+(ru(i+1)*p (i+1)/r(i+1))...
		   -(rE(i)*  ru (i)  /r(i))-  (ru(i)  *p (i)  /r(i)));
end

for i=1:nx,ph(i)=(gg-1)*(rEh(i)-0.5*(ruh(i)*ruh (i)/rh(i)));end

for i=2:nx-1					% correction step
 r(i)=r(i)-(dt/h)*(ruh(i)-ruh(i-1));
 ru(i)=ru(i)-(dt/h)*((ruh (i)^2/rh(i))-(ruh (i-1)^2/rh(i-1))+ph(i)-ph(i-1));
 rE(i)=rE(i)-(dt/h)*((rEh(i)*ruh (i)/rh(i))-(rEh(i-1)*ruh (i-1)/rh(i-1))+...
	(ruh(i)*ph (i)/rh(i))-(ruh(i-1)*ph (i-1)/rh(i-1)));
end

for i=1:nx,u(i)=ru (i)/r(i);end

for i=2:nx-1					% artificial viscosity
 ru(i)=ru(i)+artvisc*(0.5*dt/h)*( (r(i)+r(i+1))*(u(i+1)-u(i))*abs(u(i+1)-u(i))...
                              -(r(i)+r(i-1))*(u(i)-u(i-1))*abs(u(i)-u(i-1)) );

 rE(i)=rE(i)+artvisc*(0.25*dt/h)*...
                     ( (u(i+1)+u(i))*(r(i)+r(i+1))*(u(i+1)-u(i))*abs(u(i+1)-u(i))...
                      -(u(i)+u(i-1))*(r(i)+r(i-1))*(u(i)-u(i-1))*abs(u(i)-u(i-1)) );
end
time=time+dt,istep
plot(p,'b','linewidth',2),title('pressure')
if(time > 0.005)break,end
pause(0.5)
end



hold on
for i=1:nx,x(i)=(i-1)/(nx-1);end
plot(x,r,'b','linewidth',2),title('density')
axis([0, 1, 0, 1.1*max(r)]);
set(gca,'fontsize',16);set(gca,'linewidt',2)
box on
% hold on
%	plot(p,'b'),axis([0 81 0 150000])
% hold off
% end
