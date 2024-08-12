function dout = rflsm_ssa(din,t,param)
% This function applies Singular Spectrum Analysis (SSA) filter to receiver
% functions.
% input:      din -- input data
%             t -- time axis
%             h -- offset axis
%             param -- struct contains the parameters for radon transform
% output:     dout -- filtered data
%
% August, 2024, Yunfeng Chen, write the function
global ishot
disp('SSA begins')
flow = param.flow;    % minimum frequency
fhigh = param.fhigh;  % maximum frequency
rank_p=param.rank_p;  % rank value
alpha=param.alpha;    % trade-off parameter
n_iter=param.n_iter;  % iteration number
xmax=param.xmax;

% Compute the Sampling Operator S
T = ones(size(din));
[nt,ntrace] = size(din);
for ix = 1:ntrace
    a = sum(din(:,ix));
    if a==0
        T(:,ix)=0;
    end
end
% Reconstruction
dt = t(2)-t(1);
tic;
dout = reconstruction(din,T,dt,flow,fhigh,rank_p,alpha,n_iter);
toc;
disp('SSA done')
%% Plot ssa results
figdir=param.figdir;
fig=figure;
set(gcf,'Position',[100 100 1000 500],'Color','w')
subplot(121)
imagesc([0,xmax],t,din);
caxis([-0.1 0.1])
colorbar
ylim([0 30])
xlabel('Distance (km)');
ylabel('Time (sec)');
set(gca,'fontsize',14)
text(-0.2,0.98,'(a)','Units','normalized','FontSize',18)

subplot(122)
imagesc([0,xmax],t,dout);
caxis([-0.1 0.1])
colorbar
colormap(seismic(3));
ylim([0 30])
xlabel('Distance (km)');
ylabel('Time (sec)');
set(gca,'fontsize',14)
text(-0.2,0.98,'(b)','Units','normalized','FontSize',18)
figname=['rf_ssa.',num2str(ishot),'.png'];
export_fig(fig,fullfile(figdir,figname));