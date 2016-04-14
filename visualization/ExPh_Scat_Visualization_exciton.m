%% This file visualizes the exciton-phonon scattering rates
clear all; clc; fig=40;
% close all;
dir='C:\Users\Amirhossein\Documents\My Dropbox\Research\Exciton\FortranCode\CNT_17_00_4001_0200_1.5_1\';
eV=1.6e-19;


%% plot CNT energy dispersion
% FileName=[dir,'electron_conduction.dat'];
% Ec_tmp=load(FileName);
% FileName=[dir,'electron_valence.dat'];
% Ev_tmp=load(FileName);
% 
% [Nu,nkc]=size(Ec_tmp);
% Nu=Nu-1;
% k_vec=Ec_tmp(1,:);
% 
% E_c=Ec_tmp(2:Nu+1,:);
% E_v=Ev_tmp(2:Nu+1,:);
% 
% fig=fig+1; figure(fig); hold on; box on;
% plot(k_vec,E_c/eV,'-','LineWidth',3);
% plot(k_vec,E_v/eV,'-','LineWidth',3);
% axis tight;

%% plot the crossing points
FileName=[dir,'TotScatRate.dat'];
ScatRate=load(FileName);

fig=fig+1; figure(fig); box on;
semilogy(ScatRate(:,1)/eV,ScatRate(:,2)*eV,'--b','LineWidth',3); hold on;
% semilogy(ScatRateN(:,1)/eV,ScatRateN(:,2),'-b','LineWidth',3);
axis tight;