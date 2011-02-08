%% cleaning
clear all
clc

%% read necessary files
vy=dlmread('vy.txt',' ');
vz=dlmread('vz.txt',' ');
phase=dlmread('phase.txt',' ');
%% draw streamlines

% specifying the starting poitns

%[sx,sy] = meshgrid(0:20:1500,0:10:50);

sy=1:52;
sz=ones(1,52);

fig=figure()
set(gcf,'PaperUnits','centimeters')
xSize = 8; ySize = 12;
%xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 1500 400])

%axis=gca()
%get(axis)
h=streamline(vz,vy,sz,sy,[0.1, 15000])
set(h,'Color','red')
figure()
imshow(phase)