%% This file visualizes the exciton-phonon scattering process
clearvars -except phk_vec ph_disp Kcm_vec Ex0_A2 Nu dk; clc; fig=30;
% close all;
dir='C:\Users\Amirhossein\Documents\My Dropbox\Research\Exciton\FortranCode\CNT_10_00_0501_0200_1.5_1\';
eV=1.6e-19;

%% plot CNT dispersion for both phonons and excitons
[nKcm,nX]=size(Ex0_A2);
nkph=numel(phk_vec);

fig=fig+1; figure(fig); hold on; box on;
for i=1:nX
    plot(Kcm_vec,Ex0_A2(:,i)/eV,'-r.','LineWidth',3,'MarkerSize',20);
    hold on;
end;
axis tight;
for i=1:6
    plot(Kcm_vec(39)-phk_vec(:)/2,Ex0_A2(39,1)/eV+ph_disp(Nu/2,:,i)/eV,'-','LineWidth',3);
end;

%% plot the crossing points
FileName=[dir,'crossmode.dat'];
crossmode=load(FileName);

plot(crossmode(:,1)*dk,crossmode(:,2)/eV,'k.','MarkerSize',20);

% axis equal; axis tight;