function [src,pos,tshift] = rflsm_create_src(dt,nfft,rayp,baz,vp,param)
% This function creates source time function required by migration.
% 
% input:      dt -- time interval
%             nfft -- number of data points
%             rayp -- ray parameter in s/km
%             baz -- back-azimuth
%             vp -- P-wave velocity at the bottom of the model
%             param -- struct contains the parameters for radon transform
% output:     src -- source time function
%             pos -- positions of point sources
%             tshift -- time lag between point sources
%             
% August, 2024, Yunfeng Chen, write the function
% get filter in F domain 
gauss = param.gauss;
ph = param.ph;
nx = param.nx;
dx = param.dx;

gaussF = gaussFilter( dt, nfft, gauss );
% phase shift in radians
shift_i = round(ph/dt); % removed +1 from here.
p = 2*pi*(1:nfft).*shift_i./(nfft);

% apply shift
Xf = gaussF.*(cos(p) - 1i .* sin(p) );

% back into time
src = real( ifft(Xf, nfft) )/cos(2*pi*shift_i/nfft);

delf = (1/dt)/nfft;
gnorm = sum(gaussF)*delf*dt;
src = real(src(:))/gnorm;

% calcualte the positions of sources
len = (nx-1)*dx;
pos{1} = linspace(0,len,nx);

% calcualte the amount of time shift
rayp=rayp*rad2km(1);
theta= -asind(rayp*vp/(6371-200)); % incidence angle

% This part of the code is kind of ad-hoc and is used to determine the
% incidence direction. You have to set the azimuth range properly based on
% your own dataset.
if baz<=270 && baz>=90
    theta=-theta;
end

tmax=tand(theta)*len/mean(vp);
if tmax<0
    tshift=linspace(abs(tmax),0,nx);
else
    tshift=linspace(0,abs(tmax),nx);
end


