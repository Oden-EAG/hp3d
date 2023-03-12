




function plot_geom(nrcurves,Curves,nrquad,Surfaces,nrhexas,Volume,Points)


figure(1);
for i = 1:nrcurves
    id1 = Curves(i,1);
    id2 = Curves(i,2);
    v1 = Points(id1,:);
    v2 = Points(id2,:);
    v =  [v1;v2];
    plot3(v(:,1), v(:,2), v(:,3), '-sk','LineWidth',3);
    xlim([-2,2]);
    ylim([-2,2]);
    zlim([0,4]);
    hold on
end


% figure(2);
% for i = 1:nrquad
%     id1 = Surfaces(i,1);
%     id2 = Surfaces(i,2);
%     id3 = Surfaces(i,3);
%     id4 = Surfaces(i,4);
%     v1 = Points(id1,:);
%     v2 = Points(id2,:);
%     v3 = Points(id3,:);
%     v4 = Points(id4,:);
%     v =  [v1;v2;v3;v4;v1];
%     plot3(v(:,1), v(:,2), v(:,3), '-sk','LineWidth',3);
%     xlim([-2,2]);
%     ylim([-2,2]);
%     zlim([0,4]);
% end

figure;

for i=1:nrhexas
    id = Volume(i,:);
    vert=[];
    for j = 1:8
        k = id(j);
        xyz = Points(k,:);
        vert = [vert ; xyz];
    end
    faces = [1 2 3 4; 5 6 7 8; 2 3 7 6; 1 4 8 5; 1 2 6 5  ; 3 4 8 7  ];
    patch('Vertices', vert, 'Faces', faces, 'FaceColor', 'g', 'FaceAlpha',0.5, 'EdgeColor','k', 'EdgeAlpha',1,'LineWidth',3)
    hold on;
    xmin = min(Points(:,1)); xmax = max(Points(:,1));
    ymin = min(Points(:,2)); ymax = max(Points(:,2));
    zmin = min(Points(:,3)); zmax = max(Points(:,3));
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    zlim([zmin,zmax]);
    az = -30;
    el = 19;
    view(az,el)
    grid on   
end

% figure;
% 
% for i = 57:nrhexas
%     id = Volume(i,:);
%     vert=[];
%     for j = 1:8
%         k = id(j);
%         xyz = Points(k,:);
%         vert = [vert ; xyz];
%     end
%     
%     vertb = vert(1:4,:) ;
%     vertb = [vertb; vertb(1,:)];
%     vertt = vert(5:8,:) ;
%     vertt = [vertt; vertt(1,:)];
%     for j = 1:4
%         v1 = vertb(j,:);
%         v2 = vertb(j+1,:);
%         v = [v1; v2];
%         plot3(v(:,1), v(:,2), v(:,3), '-','LineWidth',3);
%         xlim([0,4]);
%         ylim([0,4]);
%         zlim([0,4]);
%         az = -30;
%         el = 19;
%         view(az,el)
%         hold on;
%     end    
%     for j = 1:4
%         v1 = vertt(j,:);
%         v2 = vertt(j+1,:);
%         v = [v1; v2];
%         plot3(v(:,1), v(:,2), v(:,3), '-','LineWidth',3);
%         hold on;
%     end    
% end
% 
