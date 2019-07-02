
clear; clc; close all

az = 54;
el = 24;

b=textread('matrixMesh.txt','%f');
NRELES = b(end);
a = b(1:end-NRELES-1);
% a=reshape(a,24,length(a)/24)';
% tic
% x=a(:,3*[0:7]+1); y=a(:,3*[0:7]+2); z=a(:,3*[0:7]+3);
% 
% figure(1)
% 
% for i=1:size(a,1)
%     
%     
% 
%     a='g';
%     plot3(x(i,[1 2]),y(i,[1 2]),z(i,[1 2]),a  ); hold on
%     plot3(x(i,[2 3]),y(i,[2 3]),z(i,[2 3]),a  ); hold on
%     plot3(x(i,[3 4]),y(i,[3 4]),z(i,[3 4]),a  ); hold on
%     plot3(x(i,[4 1]),y(i,[4 1]),z(i,[4 1]),a ); hold on
%     
%     
%     plot3(x(i,[5 6]),y(i,[5 6]),z(i,[5 6]),a  ); hold on
%     plot3(x(i,[6 7]),y(i,[6 7]),z(i,[6 7]),a ); hold on
%     plot3(x(i,[7 8]),y(i,[7 8]),z(i,[7 8]),a); hold on
%     plot3(x(i,[8 5]),y(i,[8 5]),z(i,[8 5]),a ); hold on
%     plot3(x(i,[1 5]),y(i,[1 5]),z(i,[1 5]),a  ); hold on
%     plot3(x(i,[2 6]),y(i,[2 6]),z(i,[2 6]),a ); hold on
%     plot3(x(i,[3 7]),y(i,[3 7]),z(i,[3 7]),a); hold on
%     plot3(x(i,[4 8]),y(i,[4 8]),z(i,[4 8]),a  ); hold on
%    
%     
% end

% xlabel('x');ylabel('y');zlabel('z');
% % grid on;
% grid off; axis off
% axis equal
% xlim([0 max(max(x))]);ylim([0 max(max(y))]);zlim([0 max(max(z))]);
% view(az,el)
% toc
% tic

a = b(1:end-NRELES-1);
orders = b(end-NRELES:end);
figure(2)
for i = 1:NRELES
    nord= orders(i);
    switch(nord)
        case 1
            col = 'b';
        case 2
            col = 'c';
        case 3
            col = [0 0.5 0];
        case 4
            col = 'y';
        case 5
            %orange
            col = [1 .5 0];
        case 6
            col = 'r';
        case 7
            col = 'm';
        case 8
            col = [1 0.4 0.6];
        otherwise
            disp('Unknown order') ; pause;
    end
    for j = 1:8
        x(j) = a((i-1)*24 + 3*(j-1) + 1);
        y(j) = a((i-1)*24 + 3*(j-1) + 2);
        z(j) = a((i-1)*24 + 3*(j-1) + 3);
    end
%     if (max(x)+max(y)+max(z) < 1) 
    vert  = [x(1) y(1) z(1); x(2) y(2) z(2); x(3) y(3) z(3); x(4) y(4) z(4); ...
             x(5) y(5) z(5); x(6) y(6) z(6); x(7) y(7) z(7); x(8) y(8) z(8);];
%     
    faces = [1 2 3 4; 5 6 7 8; 2 3 7 6; 1 4 8 5; 1 2 6 5  ; 3 4 8 7  ];
    patch('Vertices', vert, 'Faces', faces, 'FaceColor', col, 'FaceAlpha',0.05, 'EdgeColor',col, 'EdgeAlpha',0.5)
    hold on;    
%     end
end
xlabel('x');ylabel('y');zlabel('z');
grid off; axis off
axis equal


view(az,el)
 toc