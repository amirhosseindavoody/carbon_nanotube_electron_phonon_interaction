%% This file visualizes the exciton-phonon scattering rates
clear all; clc; fig=40;
close all;
dir='C:\Users\Amirhossein\Documents\My Dropbox\Research\Exciton\FortranCode\test\';
eV=1.6e-19;


%% plot CNT energy dispersion
FileName=[dir,'electron_conduction.dat'];
Ec_tmp=load(FileName);
FileName=[dir,'electron_valence.dat'];
Ev_tmp=load(FileName);

[Nu,nkc]=size(Ec_tmp);
Nu=Nu-1;
k_vec=Ec_tmp(1,:);

E_c=Ec_tmp(2:Nu+1,:);
E_v=Ev_tmp(2:Nu+1,:);

fig=fig+1; figure(fig); hold on; box on;
% plot(k_vec,E_c/eV,'-','LineWidth',3);
% plot(k_vec,E_v/eV,'-','LineWidth',3);
plot([-(nkc-1)/2:(nkc-1)/2],E_c/eV,'-','LineWidth',3);
plot([-(nkc-1)/2:(nkc-1)/2],E_v/eV,'-','LineWidth',3);
axis tight;

%% plot the crossing points
FileName=[dir,'ScatRateP.dat'];
ScatRateP=load(FileName);

FileName=[dir,'ScatRateN.dat'];
ScatRateN=load(FileName);

fig=fig+1; figure(fig); box on;
semilogy(ScatRateP(:,1)/eV,ScatRateP(:,2),'-r','LineWidth',3); hold on;
semilogy(ScatRateN(:,1)/eV,ScatRateN(:,2),'-b','LineWidth',3);
axis tight;