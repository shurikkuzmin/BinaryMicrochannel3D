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
     %sz=15*ones(1,2*size(vz,1)-2)/size(vz,2);
     if counter==1
          sz=6000*ones(1,2*size(vz,1)-2)/size(vz,2);
     else
         sz=18000*ones(1,2*size(vz,1)-2)/size(vz,2); 
         sz2=2000*ones(1,2*size(vz,1)-2)/size(vz,2);
     end

     %figure()
     %plot(vz(1,:))
     %figure()
     %plot(vz(25,:))
     
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
     %h=streamline(x,y,vz-0.05,vy,sz,sy,[0.1, 15000])
     if counter==1
        h=streamline(x,y,vz-0.0105,vy,sz,sy,[0.1, 15000])    
        %%Calculation of vorticity
        extract=444
        vz_ext=vz(:,extract)-0.0105
        vy_ext=vy(:,extract)
        vy_ext_plus=vy(:,extract+1)
        vy_ext_minus=vy(:,extract-1)
        
        
        for k=2:size(vz_ext)-1
            vorticity(k)=0.0*(vz(k+1)-vz(k-1))/2.0-(vy_ext_plus(k)-vy_ext_minus(k))/2.0
        end
     else
        h=streamline(x,y,vz-0.0258,vy,sz,sy,[0.1, 5000])
        h2=streamline(x,y,vz-0.0258,vy,sz2,sy,[0.1,5000])
     
        extract=1189
        vz_ext=vz(:,extract)-0.0258
        vy_ext=vy(:,extract)
        vy_ext_plus=vy(:,extract+1)
        vy_ext_minus=vy(:,extract-1)
        
        
        for k=2:size(vz_ext)-1
            vorticity(k)=0.0*(vz(k+1)-vz(k-1))/2.0-(vy_ext_plus(k)-vy_ext_minus(k))/2.0
        end
     end
     set(h,'Color','red')
     hold on
     v=0.00001
     cont=contour(x,y,phase,v,'LineWidth',2,'Color','black')
     %hold on
     xlabel('X') 
     ylabel('Y')

     figure()
     plot(abs(vorticity))
     figure()
     plot(vy_ext,'Color','green')
     figure()
     if counter==1
         plot(phase(:,1400),'Color','magenta')
     else
         plot(phase(:,750),'Color','magenta')
     end
         %imshow(phase)


end

