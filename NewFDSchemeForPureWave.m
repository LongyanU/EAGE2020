%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                May 3             %%%%%%%%%%%%%%%%%
%%%%%%%%%       Virieux 1986 Geophysics    %%%%%%%%%%%%%%%%%
%%%%%%%%% New finite difference scheme    %%%%%%%%%%%%%%%%%

clear
clc %%%%%%%
close all
% Elapsed time is 6.198007 seconds.  May 2,2017
nt=700;    % number of time steps
eps=.6;     % stability
isnap=600;    % snapshot sampling

nx=350;
nz=350;

vp=ones(nz,nx)*3000;
vp(1:nz/2,:)=2600;

vs=vp./sqrt(3);
dx=20;  % calculate space increment

x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^12*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
xs=nz/2;
zs=nz/2-40;

h=dx;

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics



p=zeros([nz nx]); Vx=p; Vz=p;
Txxx=p;
Txzz=p;
Tzzz=p;
Txzx=p;

Txx=p;
Txz=p;
Tzz=p;


Vxx=p;
Vzz=p;
Vxz=p;
Vzx=p;

Vxp=p;
Vxs=p;
Vzp=p;
Vzs=p;


coeff=[ 1.58542, -0.300275, 0.0962557, -0.0349762, 0.0121838, -0.00349652, 0.000614601];
coeff=[ 1.49255, -0.23563, 0.0624287, -0.0213218, 0.00812705, -0.00275961, 0.000669177];


tic
for it=1:nt-2,
    
    % Txx/x
    Txxx=Txx-circshift(Txx,[0 1]);
    
    
    % Txx/z
    Txxz=Txx-circshift(Txx,[1]);
    
    
    % Txz/z
    Txzz=coeff(1)*(circshift(Txz,[-1 0])-Txz)+...
        coeff(2)*(circshift(Txz,[-2 0])-circshift(Txz,[1 0]))+...
        coeff(3)*(circshift(Txz,[-3 0])-circshift(Txz,[2 0]))+...
        coeff(4)*(circshift(Txz,[-4 0])-circshift(Txz,[3 0]))+...
        coeff(5)*(circshift(Txz,[-5 0])-circshift(Txz,[4 0]))+...
        coeff(6)*(circshift(Txz,[-6 0])-circshift(Txz,[5 0]))+...
        coeff(7)*(circshift(Txz,[-7 0])-circshift(Txz,[6 0]));
    
    % Tzz/z
    Tzzz=Tzz-circshift(Tzz,[1 0]);
    
    
    % Tzz/x
    Tzzx=Tzz-circshift(Tzz,[0 1]);
    
    
    % Txz/x
    Txzx=coeff(1)*(circshift(Txz,[0 -1])-Txz)+...
        coeff(2)*(circshift(Txz,[0 -2])-circshift(Txz,[0 1]))+...
        coeff(3)*(circshift(Txz,[0 -3])-circshift(Txz,[0 2]))+...
        coeff(4)*(circshift(Txz,[0 -4])-circshift(Txz,[0 3]))+...
        coeff(5)*(circshift(Txz,[0 -5])-circshift(Txz,[0 4]))+...
        coeff(6)*(circshift(Txz,[0 -6])-circshift(Txz,[0 5]))+...
        coeff(7)*(circshift(Txz,[0 -7])-circshift(Txz,[0 6]));
    
    
    Vxp=Vxp+dt*vp.^2./(2*vp.^2-2*vs.^2).*(Txxx+Tzzx)./h;
    Vzp=Vzp+dt*vp.^2./(2*vp.^2-2*vs.^2).*(Txxz+Tzzz)./h;
    
    temp=Txzz-1./(2*vp.^2-2*vs.^2).*(vp.^2.*Tzzx-(vp.^2-2*vs.^2).*Txxx);
    Vxs=Vxs+dt*temp/h;
    
    temp=Txzx-1./(2*vp.^2-2*vs.^2).*(vp.^2.*Txxz-(vp.^2-2*vs.^2).*Tzzz);
    Vzs=Vzs+dt*temp/h;
    
    Vx=Vxp+Vxs;
    Vz=Vzp+Vzs;
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    % Vz(zs,xs)=Vz(zs,xs)+src(it);
    
    Vxx=coeff(1)*(circshift(Vx,[0 -1])-circshift(Vx,[0 0]))+...
        coeff(2)*(circshift(Vx,[0 -2])-circshift(Vx,[0 1]))+...
        coeff(3)*(circshift(Vx,[0 -3])-circshift(Vx,[0 2]))+...
        coeff(4)*(circshift(Vx,[0 -4])-circshift(Vx,[0 3]))+...
        coeff(5)*(circshift(Vx,[0 -5])-circshift(Vx,[0 4]))+...
        coeff(6)*(circshift(Vx,[0 -6])-circshift(Vx,[0 5]))+...
        coeff(7)*(circshift(Vx,[0 -7])-circshift(Vx,[0 6]));
    
    Vxz=(circshift(Vx,[0])-circshift(Vx,[1]));
    
    
    Vzz=coeff(1)*(circshift(Vz,[-1])-circshift(Vz,[0]))+...
        coeff(2)*(circshift(Vz,[-2])-circshift(Vz,[1]))+...
        coeff(3)*(circshift(Vz,[-3])-circshift(Vz,[2]))+...
        coeff(4)*(circshift(Vz,[-4])-circshift(Vz,[3]))+...
        coeff(5)*(circshift(Vz,[-5])-circshift(Vz,[4]))+...
        coeff(6)*(circshift(Vz,[-6])-circshift(Vz,[5]))+...
        coeff(7)*(circshift(Vz,[-7])-circshift(Vz,[6]));
    
    
    Vzx=(circshift(Vz,[0 0])-circshift(Vz,[0 1]));
    
    
    
    Txx=Txx+dt*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*(vs.^2).*(Vxz+Vzx)/h;
    
    Txx(zs,xs)=Txx(zs,xs)+src(it);
    Tzz(zs,xs)=Tzz(zs,xs)+src(it);
    
    
    
    
%     if rem(it,isnap)== 0,
%         imagesc(x,z,Vx), axis equal
%         colormap gray
%         xlabel('x'),ylabel('z')
%         title(sprintf(' Time step: % i - Max ampl: % g ',it,max(max(Vx))))
%         drawnow
%     end
    
%     seis_record(it,:)=p(60,:);
%     if it==500
%         pp1=p;
%     elseif it==1000
%         pp2=p;
%     elseif it==1500
%         pp3=p;
%     elseif it==2000
%         pp4=p;
%     elseif it==2500
%         pp5=p;
%     end
    
end


toc
%
save('NewFDSchemeForPureWave.mat')

figure;imagesc(Vzp(45:end-45,45:end-45)+Vxp(45:end-45,45:end-45),[-6*10^2 3*10^2])
ylabel('z/dz')
xlabel('x/dx')
title('')

figure;imagesc(Vzs(45:end-45,45:end-45)+Vxs(45:end-45,45:end-45),[-6*10^2 3*10^2])
title('')
ylabel('z/dz')
xlabel('x/dx')

figure;imagesc(Vx(45:end-45,45:end-45),[-6*10^2 3*10^2])
ylabel('z/dz')
xlabel('x/dx')
title('')

figure;imagesc(Vz(45:end-45,45:end-45),[-6*10^2 3*10^2])
ylabel('z/dz')
xlabel('x/dx')
title('')