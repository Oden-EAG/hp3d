
clear; clc; close all
%
for j = 0:-1:0
    tic
    imesh = 1*j;
    filename=strcat('mesh_',num2str(imesh),'.txt');
    b=textread(filename,'%f');
    NRELES = b(end);
    a = b(1:end-NRELES-1);
    orders = b(end-NRELES:end);
    
    
    fprintf(' Number of elements = %d \n ',NRELES);
    
    
    %plane cut
    %normal vector
    % n = [0.18 0.59 -0.78];
    n = [0.5 0.5 -1];
    origin = [0.5 0.5 0.5]';
    
    d = n*origin;
    
    
    Fig1 = figure;
    set(Fig1, 'Position', [100, 100,1024,776]);
    for i = 1:NRELES
        nord= orders(i);
        %     fprintf('Element %d has order %d \n ',i,nord);
        switch(nord)
            case 1
                col = 'c';
            case 2
                col = 'b';
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
        if (n(1)*max(x)+n(2)*max(y)+n(3)*max(z) >= d)
            vert  = [x(1) y(1) z(1); x(2) y(2) z(2); x(3) y(3) z(3); x(4) y(4) z(4); ...
                     x(5) y(5) z(5); x(6) y(6) z(6); x(7) y(7) z(7); x(8) y(8) z(8);];
            %
            faces = [1 2 3 4; 5 6 7 8; 2 3 7 6; 1 4 8 5; 1 2 6 5  ; 3 4 8 7  ];
            patch('Vertices', vert, 'Faces', faces, 'FaceColor', col, 'FaceAlpha',0.8, 'EdgeColor','k', 'EdgeAlpha',1)
            hold on;
        end
    end
    
    h(1) = plot(NaN,NaN,'s','MarkerFaceColor','c','MarkerSize',20,'LineWidth',2);
    h(2) = plot(NaN,NaN,'s','MarkerFaceColor','b','MarkerSize',20,'LineWidth',2);
    h(3) = plot(NaN,NaN,'s','MarkerFaceColor','[0 0.5 0]','MarkerSize',20,'LineWidth',2);
    h(4) = plot(NaN,NaN,'s','MarkerFaceColor','y','MarkerSize',20,'LineWidth',2);
    h(5) = plot(NaN,NaN,'s','MarkerFaceColor','[1 .5 0]','MarkerSize',20,'LineWidth',2);
    h(6) = plot(NaN,NaN,'s','MarkerFaceColor','r','MarkerSize',20,'LineWidth',2);
    h(7) = plot(NaN,NaN,'s','MarkerFaceColor','m','MarkerSize',20,'LineWidth',2);
    h(8) = plot(NaN,NaN,'s','MarkerFaceColor','[1 0.4 0.6]','MarkerSize',20,'LineWidth',2);
    
    hh = legend(h, 'p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','location','southoutside');
    set(hh,'Orientation','horizontal')
    set(hh,'FontSize',20);
    
    
    xlabel('x');ylabel('y');zlabel('z');
    grid off; axis off
    axis equal
    az = -30;
    el = 19;
    view(az,el)
    set(gca,'FontSize',20)
    set(gca,'position',[0.01 0.01 0.98 0.98],'units','normalized')
    file1=strcat('primal_mesh',num2str(imesh),'.png');
    set(gcf,'PaperPositionMode','auto')
    print(Fig1,file1,'-dpng','-r0');    
    figure(1); clf
    toc
end
