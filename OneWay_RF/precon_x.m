function m=precon_x(u)
%  Preconditioning operator, which essentially smooth the model laterally
%        |2 1 0 0 0|
%        |1 2 1 0 0|  
%  P=0.5 |0 1 2 1 0|
%        |0 0 1 2 1|
%        |0 0 0 1 2|
% Oct. 17, 2020, Yunfeng Chen, UofA
[nz,nx,ns]=size(u);
m=zeros(nz,nx,ns);
P=0.5*[1 2 1];

for iz=1:nz
    for is=1:ns
        utmp=squeeze(u(iz,:,is));
        mtmp=conv(utmp,P,'same');
        m(iz,:,is)=mtmp;
    end
end