clear; clc; close all
% -------------  read from source -------------------------
x=load(strcat(['cube.node'])); x=x(:,3:5);
conn=load(strcat(['cube_volu.elem']));
dom=conn(:,3)-1; elim = conn(:,8:13); conn=conn(:,4:7);
n=size(conn,1); lines=zeros(n*6,2); areas=zeros(n*4,3);

% max(max(conn))

% remove extra nodes and shift ordering in connectivity
elim = sort(unique(elim(:)),'descend');
x(elim,:)=[];
for i=1:length(elim)
    conn(conn > elim(i))=conn(conn > elim(i))-1;
end

for i=1:n
    lines((i-1)*6+1,:)=sort([conn(i,1) conn(i,2)]);
    lines((i-1)*6+2,:)=sort([conn(i,2) conn(i,3)]);
    lines((i-1)*6+3,:)=sort([conn(i,3) conn(i,4)]);
    lines((i-1)*6+4,:)=sort([conn(i,4) conn(i,1)]);
    lines((i-1)*6+5,:)=sort([conn(i,2) conn(i,4)]);
    lines((i-1)*6+6,:)=sort([conn(i,1) conn(i,3)]);
end
lines_b=lines; lines=unique(lines,'rows');

for i=1:n
    areas((i-1)*4+1,:)=sort([conn(i,1) conn(i,2) conn(i,3)]);
    areas((i-1)*4+2,:)=sort([conn(i,1) conn(i,2) conn(i,4)]);
    areas((i-1)*4+3,:)=sort([conn(i,1) conn(i,3) conn(i,4)]);
    areas((i-1)*4+4,:)=sort([conn(i,2) conn(i,3) conn(i,4)]);
end

areas=unique(areas,'rows');

numpoints = size(x,1)
numsegs   = size(lines,1)
numtris   = size(areas,1)
numtets   = size(conn,1)

% ------------------   write geometry file ---------------

[fid,message]= fopen('ArteryGeom.txt','w+');
fprintf(fid,'%s\n','3   3     NDIM,MANDIM');fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','0         no surfaces');fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','3         NRDOMAIN');fprintf(fid,'%s\n','');
fprintf(fid,'%s\n',strcat(num2str(size(x,1)), '          NRPOINT'));fprintf(fid,'%s\n','');

for i=1:numpoints
% for i=1:100
    fprintf(fid,'%s','Regular             point');fprintf(fid,'%s',' ');fprintf(fid,'%s\n',num2str(i));
    fprintf(fid,'%s',strcat(num2str(x(i,1),'%2.14f'),'d0')); fprintf(fid,'%s',' '); fprintf(fid,'%s',strcat(num2str(x(i,2),'%2.14f'),'d0'));
    fprintf(fid,'%s',' ');fprintf(fid,'%s\n',strcat(num2str(x(i,3),'%2.14f'),'d0'));
    fprintf(fid,'%s\n','');
end

fprintf(fid,'%s\n',strcat(num2str(size(lines,1)), '          NRCURVE'));fprintf(fid,'%s\n','');

for i=1:numsegs
% for i=1:100
    fprintf(fid,'%s','Seglin              curve');fprintf(fid,'%s',' ');fprintf(fid,'%s\n',num2str(i));
    fprintf(fid,'%s',num2str(lines(i,1)));      fprintf(fid,'%s',' ');     fprintf(fid,'%s\n',num2str(lines(i,2)));
    fprintf(fid,'%s\n','');
end

% fprintf(fid,'%s\n','0             NRTRIAN');fprintf(fid,'%s\n','');
fprintf(fid,'%s\n',strcat(num2str(size(areas,1)), '          NRTRIAN'));fprintf(fid,'%s\n','');

for i=1:numtris
% for i=1:100
    fprintf(fid,'%s','PlaneTri            triangle');fprintf(fid,'%s',' ');fprintf(fid,'%s\n',num2str(i));
    fprintf(fid,'%s',num2str(areas(i,1))); fprintf(fid,'%s',' '); fprintf(fid,'%s',num2str(areas(i,2)));fprintf(fid,'%s',' ');
    fprintf(fid,'%s\n',num2str(areas(i,3)));
    fprintf(fid,'%s\n','');
end

fprintf(fid,'%s\n','0             NRRECTA');fprintf(fid,'%s\n','');
% fprintf(fid,'%s\n',strcat(num2str(size(areas,1)), '          NRRECTA'));fprintf(fid,'%s\n','');

fprintf(fid,'%s\n','0     NRPRISM');fprintf(fid,'%s\n','');
% fprintf(fid,'%s\n',strcat(num2str(size(conn,1)), '          NRHEXAS'));fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','0     NRHEXAS');fprintf(fid,'%s\n','');


% for i=1:size(conn,1)
%     fprintf(fid,'%s','Linear              hexa');fprintf(fid,'%s',' ');fprintf(fid,'%s\n',num2str(i));
%     fprintf(fid,'%s',num2str(i));fprintf(fid,'%s','    ');
%     fprintf(fid,'%s',num2str(conn(i,1))); fprintf(fid,'%s',' '); fprintf(fid,'%s',num2str(conn(i,2)));fprintf(fid,'%s',' ');
%     fprintf(fid,'%s',num2str(conn(i,3))); fprintf(fid,'%s',' '); fprintf(fid,'%s',num2str(conn(i,4)));fprintf(fid,'%s',' ');
%     fprintf(fid,'%s',num2str(conn(i,5))); fprintf(fid,'%s',' '); fprintf(fid,'%s',num2str(conn(i,6)));fprintf(fid,'%s',' ');
%     fprintf(fid,'%s',num2str(conn(i,7))); fprintf(fid,'%s',' '); fprintf(fid,'%s\n',num2str(conn(i,8)));
%     fprintf(fid,'%s\n','');
% end

% fprintf(fid,'%s\n','0     NRTETRA');fprintf(fid,'%s\n','');
fprintf(fid,'%s\n',strcat(num2str(size(conn,1)), '          NRTETRA'));fprintf(fid,'%s\n','');

for i=1:numtets
% for i=1:100
    fprintf(fid,'%s','Linear              tetr');fprintf(fid,'%s',' ');fprintf(fid,'%s\n',num2str(i));
    fprintf(fid,'%s',num2str(dom(i)));fprintf(fid,'%s','    ');
    fprintf(fid,'%s',num2str(conn(i,1))); fprintf(fid,'%s',' '); fprintf(fid,'%s',num2str(conn(i,2)));fprintf(fid,'%s',' ');
    fprintf(fid,'%s',num2str(conn(i,3))); fprintf(fid,'%s',' '); fprintf(fid,'%s',num2str(conn(i,4)));fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n','');
end

fprintf(fid,'%s\n','0     NRPYRAM');fprintf(fid,'%s\n','');

fclose('all');
% 
% % -------------------  plot  elem  --------------------------------
% 
% % for i=1:size(conn,1)   
% for i=1:5000 
%     
%     for j=1:6
%         tx=[x(lines_b((i-1)*6+j,1),1) x(lines_b((i-1)*6+j,2),1)];
%         ty=[x(lines_b((i-1)*6+j,1),2) x(lines_b((i-1)*6+j,2),2)];
%         tz=[x(lines_b((i-1)*6+j,1),3) x(lines_b((i-1)*6+j,2),3)];
%         plot3(tx([1 2]),ty([1 2]),tz([1 2])  ); hold on
%     end
%     
% %     tx=[x(conn(i,1),1) x(conn(i,2),1)];
% %     ty=[x(conn(i,1),2) x(conn(i,2),2)];
% %     tz=[x(conn(i,1),3) x(conn(i,2),3)];
% %     plot3(tx([1 2]),ty([1 2]),tz([1 2]),'k:','linewidth',1.5  ); hold on
% %     
% %     tx=[x(conn(i,2),1) x(conn(i,3),1)];
% %     ty=[x(conn(i,2),2) x(conn(i,3),2)];
% %     tz=[x(conn(i,2),3) x(conn(i,3),3)];
% %     plot3(tx([1 2]),ty([1 2]),tz([1 2]),'k:','linewidth',1.5  ); hold on
% %     
% %     tx=[x(conn(i,3),1) x(conn(i,4),1)];
% %     ty=[x(conn(i,3),2) x(conn(i,4),2)];
% %     tz=[x(conn(i,3),3) x(conn(i,4),3)];
% %     plot3(tx([1 2]),ty([1 2]),tz([1 2]),'k:','linewidth',1.5  ); hold on
% %     
% %     tx=[x(conn(i,4),1) x(conn(i,1),1)];
% %     ty=[x(conn(i,4),2) x(conn(i,1),2)];
% %     tz=[x(conn(i,4),3) x(conn(i,1),3)];
% %     plot3(tx([1 2]),ty([1 2]),tz([1 2]),'k:','linewidth',1.5  ); hold on
%     
% end
% xlabel('x');ylabel('y');zlabel('z');
% % view([120 50]); grid on; axis equal
% % xlim([0 2]);ylim([0 1]);zlim([0 1]);

%------------------------  plot BC faces   -----------------------------------------------------
% bc=0;j=0;
% for i=1:size(areas,1)
% 
%     tx=[x(areas(i,1),1) x(areas(i,2),1) x(areas(i,3),1) x(areas(i,4),1)];
%     ty=[x(areas(i,1),2) x(areas(i,2),2) x(areas(i,3),2) x(areas(i,4),2)];
%     tz=[x(areas(i,1),3) x(areas(i,2),3) x(areas(i,3),3) x(areas(i,4),3)];
% 
%     if (    (tx(1)==0 && tx(2)==0  && tx(3)==0 && tx(4)==0)) % back
%  %    if (    (tx(1)==4 && tx(2)==4  && tx(3)==4 && tx(4)==4)) % front
%  %      if (    (ty(1)==0 && ty(2)==0  && ty(3)==0 && ty(4)==0)) % left
% %         if (    (ty(1)==1 && ty(2)==1  && ty(3)==1 && ty(4)==1)) % right
% %        if (    (tz(1)==0 && tz(2)==0  && tz(3)==0 && tz(4)==0)) % bottom
% %     if (    (tz(1)==1 && tz(2)==1  && tz(3)==1 && tz(4)==1)) % top
% 
%         j=j+1;      bc(j)=i;
%         plot3(tx([1 2]),ty([1 2]),tz([1 2])  ); hold on
%         plot3(tx([2 3]),ty([2 3]),tz([2 3]) );
%         plot3(tx([3 4]),ty([3 4]),tz([3 4]) );
%         plot3(tx([4 1]),ty([4 1]),tz([4 1]) );
%     end
% 
% end
% xlabel('x');ylabel('y');zlabel('z');
% view([120 50]); grid on; axis equal
% xlim([0 2]);ylim([0 1]);zlim([0 1]);
% bc

% --------------- write BC matrix ------------------------------------------------

% BC=zeros(size(conn,1),6);  f=zeros(6,4);
% for i=1:size(conn,1)
% 
%     f(1,:)=[conn(i,1) conn(i,2) conn(i,3) conn(i,4)];
%     f(2,:)=[conn(i,5) conn(i,6) conn(i,7) conn(i,8)];
%     f(3,:)=[conn(i,5) conn(i,6) conn(i,2) conn(i,1)];
%     f(4,:)=[conn(i,6) conn(i,2) conn(i,3) conn(i,7)];
%     f(5,:)=[conn(i,8) conn(i,7) conn(i,3) conn(i,4)];
%     f(6,:)=[conn(i,5) conn(i,1) conn(i,4) conn(i,8)];
% 
%     for j=1:6
% 
%         %%% all DBC
% %         if (abs(norm( x(f(j,:),1)' - [0 0 0 0] ))<1e-10 || ... % back
% %                 abs(norm( x(f(j,:),1)' - [2 2 2 2] ))<1e-10 || ... % front
% %                 abs(norm( x(f(j,:),2)' - [0 0 0 0] ))<1e-10 || ... % left
% %                 abs(norm( x(f(j,:),2)' - [1 1 1 1] ))<1e-10 || ... % right
% %                 abs(norm( x(f(j,:),3)' - [0 0 0 0] ))<1e-10 || ... % bottom
% %                 abs(norm( x(f(j,:),3)' - [1 1 1 1] ))<1e-10)       % top
% %             
%             
% %         %%% front/back
% %         if (abs(norm( x(f(j,:),1)' - [0 0 0 0] ))<1e-10 || ... % back
% %                 abs(norm( x(f(j,:),1)' - [2 2 2 2] ))<1e-10)   % front
%             
%             %         %%% /back
%         if (abs(norm( x(f(j,:),1)' - [0 0 0 0] ))<1e-10 ) % back
% 
%             BC(i,j)=1;
%         end
%     end
% end
% BC=BC'
