%% cleaning
clear all
clc

%% read necessary files
vy=dlmread('vely.dat',' ');
vz=dlmread('velz.dat',' ');
phase=dlmread('phase.dat',' ');
%% draw streamlines

% specifying the starting poitns

%[sx,sy] = meshgrid(0:20:1500,0:10:50);

sy1=1.0/(size(vz,1))*(2:size(vz,1)-1);
sz1=15*ones(1,size(vz,1)-2)/size(vz,2);

sy2=1.0/(size(vz,1))*(2:size(vz,1)-1);
sz2=11*ones(1,size(vz,1)-2);
%sy=1:size(vz,1)
%sz=ones(1,size(vz,1))

%sy=[0.1 0.2]
%sz=[1.0 1.0]

fig=figure()

%subplot(211)
set(gcf,'PaperUnits','centimeters')
xSize = 8; ySize = 12;
set(gcf,'Position',[0 0 1500 400])
[x,y]=meshgrid(1.0*15/size(vz,2)*(1:1500),1.0/size(vz,1)*(1:52));
h=streamline(x,y,vz,vy,sz1,sy1,[0.1 15000]);
hold on
v=0.00001
cont=contour(x,y,phase,v,'LineWidth',2,'Color','black')
%hold on
xlabel('X')
ylabel('Y')
%subplot(212)
%axis([0,15,0,1])
%xscale([0 15])
%set(l,'Linewidth',3)
%quiver(vz,vy)
%set(axis,'XLabel','Coorx')
%set(axis,'YLabel','Coory')
%h=streamline(x,y,vz,vy,sz2,sy2,[0.1 15000])
%xlim([0 15])
%set(axis,'DataAspectRatio',[10 1 1])
%set(h,'Color','red')
%axes('Position',[0 15 0 1])
%figure()
%imshow(vy)
legend('Streamlines','Phase contour')
print -depsc2 -tiff streamlines.eps
%saveas(fig,'streamlines.eps',format='eps')