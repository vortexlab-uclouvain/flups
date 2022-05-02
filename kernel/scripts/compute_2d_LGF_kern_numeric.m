%{
Compute Green lattice function 3D
%}

clear all; close all; clc;

%size of the kernel
N = 32;
G = zeros(N,N);
filename = ['LGF_2d_sym_acc12_' num2str(N)];

%% close points computed using integral formulation

i = sqrt(-1);
tic;
for n2 = 0:N-1
    for n1 = 0:N-1
        disp([num2str(n1) ' ' num2str(n2)]);
        fun = @(t1,t2) (cos(t1*n1+t2*n2)-1.0)./(4.0*(sin(t1/2.0).^2) + 4.0*(sin(t2/2.0).^2.0));
        [I1 err1] = quad2d(fun,-pi,0,-pi,0,'RelTol',1e-12,'AbsTol',0,'MaxFunEvals',99000,'Singular',true);
        [I2 err2] = quad2d(fun, 0,pi,-pi,0,'RelTol',1e-12,'AbsTol',0,'MaxFunEvals',99000,'Singular',true);
        [I3 err3] = quad2d(fun,-pi,0,0,pi,'RelTol',1e-12,'AbsTol',0,'MaxFunEvals',99000,'Singular',true);
        [I4 err4] = quad2d(fun, 0,pi,0,pi,'RelTol',1e-12,'AbsTol',0,'MaxFunEvals',99000,'Singular',true);
        disp(['estimated error = ' num2str((err1+err2+err3+err4)./(4*pi*pi))]);
        G(n1+1,n2+1) = 1./(4*pi*pi).*(I1+I2+I3+I4);
        toc
    end
end
toc

% save it, just in case
save([filename 'mat'],'G');

%% write that shit
exist(['./' filename '.ker'], 'file')
if (exist(['./' filename '.ker'], 'file') == 2)
    delete(filename);
end
% save everything
fileID = fopen([filename '.ker'],'w');
for n2 = 0:N-1
    for n1 = 0:N-1
        fwrite(fileID,G(n1+1,n2+1),'double');
    end
end
fclose(fileID)

%% figure, just for kidding
figure;
surf(G(:,:,1),'edgecolor','none')

figure; hold all;
plot(squeeze(G(:,1,1)),'.-');