function [mig,dp] = rflsm_lsm(d,t,vel,vel_s,src,pos,tshift,param)
% This function applies least-squares migration to receiver functions
% 
% input:      d -- time-shifted receiver functions
%             t -- time axis
%             vel -- P velocity model
%             vel_s -- S velocity model
%             src -- source time function
%             pos -- positions of point sources
%             tshift -- time lag between point sources
%             param -- struct contains the parameters for radon transform
% output:     mig -- migration image
%             dp -- predicted receiver functions
%             
% August, 2024, Yunfeng Chen, write the function
global ishot

disp('Preconditioned least-squares migration begins')

% parameters for forward propagation
fpeak = 10.; %not used
dt = t(2)-t(1);
nt = length(t);
dx = param.dx;
dz = param.dz;
x = param.x;
z = param.z;
nx = param.nx;
nz = param.nz;
xmax=param.xmax;
ph=abs(t(1));
bc = 1;
flow = param.flow;
fhigh = param.fhigh;
itermax=param.itermax;
mu=param.mu;
save_wavefield=0;

% sampling operator
S=repmat(any(d),size(d,1),1);
S=S(:);
if_cg=1;
P = @(u) precon_x(u);
Pt= P;
L = @(m) S.*kvssmfor_rf(m,vel,vel_s,pos,[],tshift-ph, ...
    nt,dt,dx,dz,flow,fhigh,bc,'P',src,save_wavefield,if_cg);
Lt= @(d) kvssmadj_rf(S.*d,vel,vel_s,pos,[],tshift-ph, ...
    nt,dt,dx,dz,flow,fhigh,bc,'P',src,save_wavefield,if_cg);
LP= @(u) L(P(u));
PtLt= @(d) Pt(Lt(d));

% dot product test
%     m1=randn(nz,nx,1);
%     d2=randn(nt,nx,1);
%     m1=m1(:);
%     d2=d2(:);
%     d1=LP(m1);
%     m2= PtLt(d2);
%     dot1=sum(sum(sum(d1.*d2)))
%     dot2=sum(sum(sum(m1.*m2)))

A = @(u) PtLt(LP(u))+mu*u;
b = PtLt(d(:));
tic;
[utmp,flag,relres,iter,resvec] = pcg(A,b,[],itermax);
toc;
utmp=reshape(utmp,nz,nx);
mtmp=P(utmp);
dp = L(mtmp(:));
mig = mtmp;
disp('Preconditioned least-squares migration done')

% plot results
figdir = param.figdir;
fig=figure;
set(gcf,'Position',[100 100 1800 600],'color','w')
subplot(131)
imagesc(x,t, reshape(d/max(d(:)),nt,nx));
xlim([0 xmax])
caxis([-1 1]);
colorbar
xlabel('Distance (km)');
ylabel('Time (sec)');
title('Time-shifted RF')
set(gca,'fontsize',14)
text(-0.2,0.98,'(a)','Units','normalized','FontSize',18)

subplot(132)
imagesc(x,z, mean(mtmp,3));hold on;
xlim([0 xmax])
colormap(seismic(3));
colorbar;
cmax=rms(abs(mtmp(:)));
caxis([-3*cmax 3*cmax]);
xlabel('Distance (km)');
ylabel('Depth (km)');
title('LSM image')
set(gca,'fontsize',14)
text(-0.2,0.98,'(b)','Units','normalized','FontSize',18)

subplot(133)
imagesc(x,t, reshape(dp/max(dp(:)),nt,nx));
xlim([0 xmax])
xlabel('Distance (km)');
ylabel('Time (sec)');
title('Predicted data')
set(gca,'fontsize',14)
caxis([-1 1]);
colorbar
text(-0.2,0.98,'(c)','Units','normalized','FontSize',18)
figname=['migls.',num2str(ishot),'.png'];
export_fig(fig,fullfile(figdir,figname));