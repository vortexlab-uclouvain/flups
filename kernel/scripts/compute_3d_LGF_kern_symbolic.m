%{
Compute Green lattice function 3D
%}

clear all; close all; clc;


%size of the kernel
N = 64;
G = zeros(N,N,N);
filename = ['LGF_3d_sym_acc12_' num2str(N)];

%% close points computed using integral formulation
tic;
for n3=0:N-1
    for n2 = 0:N-1
        for n1 = 0:N-1
            disp([num2str(n1) ' ' num2str(n2) ' ' num2str(n3)]);
            fun = @(t) besseli(n1,2*t,1).*besseli(n2,2*t,1).*besseli(n3,2*t,1);
            G(n1+1,n2+1,n3+1) = quadgk(fun,0,inf,'RelTol', 1e-12, 'AbsTol', 0);
            toc
        end
    end
end
toc

% save it, just in case
save([filename 'mat'],'G');

%% write that shit
if (exist([filename '.ker'], 'file') == 2)
    delete(filename);
end
% save everything
fileID = fopen([filename '.ker'],'w');
for n3=0:N-1
    for n2 = 0:N-1
        for n1 = 0:N-1
            fwrite(fileID,G(n1+1,n2+1,n3+1),'double');
        end
    end
end
fclose(fileID)

%% figure, just for kidding
figure;
surf(G(:,:,1),'edgecolor','none')

figure; hold all;
plot(squeeze(G(:,1,1)),'.-');