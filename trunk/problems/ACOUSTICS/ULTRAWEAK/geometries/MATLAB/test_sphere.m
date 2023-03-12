
clear all

Points = [-1 -1 -1 ; -1 -1 1 ; -1 1 -1 ; -1 1 1 ;  1 -1 -1 ; 1 -1 1 ; 1 1 -1 ; 1 1 1 ] ;

Points = [0.125*Points ; sqrt(3)/12*Points];


Curves = [

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


% 
% figure;
% for i = 1:nrhexas
%     id = Volume(i,:);
%     vert=[];
%     for j = 1:8
%         k = id(j);
%         xyz = Points(k,:);
%         vert = [vert ; xyz];
%     end
%     faces = [1 2 3 4; 5 6 7 8; 2 3 7 6; 1 4 8 5; 1 2 6 5  ; 3 4 8 7  ];
%     patch('Vertices', vert, 'Faces', faces, 'FaceColor', 'g', 'FaceAlpha',0.5, 'EdgeColor','k', 'EdgeAlpha',1,'LineWidth',3)
%     hold on;
%     xlim([0,4]);
%     ylim([0,4]);
%     zlim([0,4]);
%     az = -30;
%     el = 19;
%     view(az,el)
%     
% end