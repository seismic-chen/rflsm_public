function param = rflsm_init_param()
% This function initiates parameters used by the main function

% processing flag
param.is_radon = true;
param.is_ssa = true;

if ~param.is_radon && ~param.is_ssa
    suffix = '_no_processing';
end
if param.is_radon && ~param.is_ssa
    suffix = '_no_ssa';
end
if param.is_radon && param.is_ssa
    suffix = '';
end

param.datadir = './decon/';
figdir = ['./figures',suffix];
param.figdir = figdir;


if ~exist(figdir)
    mkdir(figdir)
end

% parameters for velocity model
param.dx=4;
param.dz=0.5;
param.xmax=551;
param.ymax=200;
param.xpad=300;
param.x=0-param.xpad:param.dx:param.xmax+param.xpad;
param.z=0:param.dz:param.ymax;
param.nx=length(param.x);
param.nz=length(param.z);

% parameters for quality control
param.flow=0.05;
param.fhigh=1.2;
param.minsnr = 0.5;
param.mintrace = 30;

% parameters for profile location
param.lon1 = 85.80837236;
param.lat1 = 29.2673;
param.lon2 = 83.87268109;
param.lat2 = 33.9673;

% parameter for receiver function deconvolution
param.gauss=2.5;
param.itmax = 400;
param.minderr = 0.001;
param.ph = 5; % phase delay
param.VB = 0; % Allow verbose output

% parameter for radon transform
param.pmax = 2;
param.dp = 0.005;
param.N1 = 15;  % CG Iterations (Internal loop);
param.N2 = 1;   % Update of weights for the sparse solution: N1 = 1 LS;  N2 > 3 for High Res (Sparse) solution

% parameter for sigular spectrum analysis
param.rank_p=8;
param.alpha=0.9;
param.n_iter=20;

% parameter for least-squares migration
param.itermax = 10;
param.mu = 0.1;
