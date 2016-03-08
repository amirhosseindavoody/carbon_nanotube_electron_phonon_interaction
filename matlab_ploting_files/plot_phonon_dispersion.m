%%
% size(a)
% size(b)
% plot(a,b(:,1))

clc
close all;
clear all;
dir='C:\Users\amirhossein\Downloads\';
eV=1.6e-19;

FileName=[dir,'phonon_dispersion_10.dat'];
tmp=load(FileName);
nc=size(tmp,1);
x = tmp(1,:);
y = tmp(2:nc,:);

% x=linspace(-10,10,47);
% y=x.^2;
figure
plot(x,y);

% figure
% plot(a,b(:,1));