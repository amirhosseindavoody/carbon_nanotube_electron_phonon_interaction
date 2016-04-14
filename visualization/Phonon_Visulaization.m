%% This file visualizes the results of the fortran program for CNT phonon dispersion
clearvars -except Kcm_vec Ex0_A2; clc; fig=20;
% close all;
dir='C:\Users\Amirhossein\Documents\My Dropbox\Research\Exciton\FortranCode\CNT_10_00_0501_0200_1.5_1\';
eV=1.6e-19;

%% plot CNT dispersion
FileName=[dir,'phonon_dispersion.dat'];
tmp=load(FileName);
[nrow,nkc]=size(tmp);
phk_vec=tmp(1,:);
dk=phk_vec(2)-phk_vec(1);
tmp1(1:nrow-1,1:nkc)=tmp(2:nrow,1:nkc);
Nu=(nrow-1)/6;

fig=fig+1; figure(fig); hold on; box on;
for i=1:6
    for mu=1:Nu
        ph_disp(mu,1:nkc,i)=tmp(1+(i-1)*Nu+mu,1:nkc);
        plot(phk_vec(1:nkc),ph_disp(mu,1:nkc,i)/eV,'-','LineWidth',2);
    end;
end;
axis tight;

for i=1:6
    for mu=Nu/2
        ph_disp(mu,1:nkc,i)=tmp(1+(i-1)*Nu+mu,1:nkc);
        plot(phk_vec(1:nkc),ph_disp(mu,1:nkc,i)/eV,'-r','LineWidth',6);
    end;
end;

% plot(phk_vec((nkc+1)/2:nkc),ph_disp(:,(nkc+1)/2:nkc)/eV,'-','LineWidth',2);
% axis tight;