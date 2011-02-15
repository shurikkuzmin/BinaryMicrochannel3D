%% cleaning
clear all
clc

%% draw streamlines

% specifying the starting poitns

%[sx,sy] = meshgrid(0:20:1500,0:10:50);

names={'Force0000002/4/','Force0000002/8/'}

for counter=1:2  
     %% read necessary files
     filenamey=strcat(names{counter},'vy.txt')
     filenamez=strcat(names{counter},'vz.txt')
     filenamephase=strcat(names{counter},'phase.txt')
     vy=dlmread(filenamey,' ');
     vz=dlmread(filenamez,' ');
     phase=dlmread(filenamephase,' ');
     size(vy)
     size(vz)
     size(phase)
     %sy=1:52;
     %sz=ones(1,52);
     
     sy=1.0/(2*size(vz,1))*(2:2*size(vz,1)-1); 
     sz=15*ones(1,2*size(vz,1)-2)/size(vz,2);

     
     fig=figure()
     set(gcf,'PaperUnits','centimeters')
     xSize = 8; ySize = 12;
     %xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
     %set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
     set(gcf,'Position',[0 0 1500 400])

     %axis=gca()
     %get(axis)
     %% Visualization part
     [x,y]=meshgrid(1.0*15/size(vz,2)*(1:1500),1.0/size(vz,1)*(1:52));
     h=streamline(x,y,vz,vy,sz,sy,[0.1, 15000])
     set(h,'Color','red')
     hold on
     v=0.00001
     cont=contour(x,y,phase,v,'LineWidth',2,'Color','black')
     %hold on
     xlabel('X') 
     ylabel('Y')

     %figure()
     %imshow(phase)


end

