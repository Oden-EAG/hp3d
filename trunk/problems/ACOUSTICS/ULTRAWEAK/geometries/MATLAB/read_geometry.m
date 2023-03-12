clear; clc; close all
format long
filename=strcat('hemisphere');
a=fileread(filename);


% read file
fid = fopen('hemisphere');
FC = textscan(fid, '%s');
a = FC{1};
fclose(fid);

%------------------------------------------------
%                READ POINTS
%------------------------------------------------
%number of points
nrpoints = a(21);
nrpoints = str2num(nrpoints{1});
b = a(23:end);
Points = [];
j=4;
for i = 1:nrpoints
    for k = 1:3
        y = b(j+k-1);
        y = y{1};
        y = y(1:end-2);
        x(k) = str2double(y);
    end
    Points(i,:) = x;
    j = j+6;
end
j = j-3;
%------------------------------------------------
%               READ CURVES
%------------------------------------------------
% straight lines 
c = b(j:end);
nrcurves = c(1);
nrcurves = str2num(nrcurves{1});
clear x
j = 5;
for i = 1:nrcurves
    for k = 1:2
        y = c(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Curves(i,:) = x;
    str = c(j+3);
    str = str{1};
    if isequal(str, 'Seglin')
        j = j+5;
    else
      break 
    end    
end
nrlines = i;
% 
% curved lines
d = c(j+3:end);
clear x
j = 3;
for i = nrlines+1:nrcurves
    for k = 1:2
        y = d(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Curves(i,:) = x;
    str = d(j+4);
    str = str{1};
    if isequal(str, '1SurfsCur')
        j = j+6;
    else
      break 
    end    
end
nrcurvlines = i-nrlines;

% Rotation    
theta(1)= 0; 
theta(2)= 0;
theta(3)=90;
rad   = theta*pi/180;
sina(1:3) = sin(rad(1:3));
cosa(1:3) = cos(rad(1:3));
    
% Rz = [cosa(3) -sina(3) 0 ; sina(3) cosa(3) 0 ; 0 0 1];
Rz = [0 -1 0 ; 1 0 0 ; 0 0 1];    
NewPoints = (Rz*Points')';

NewPoints


% Construct ONLY new curves 
for i = 1:nrcurves
    flag = 0;
    id1 = Curves(i,1);
    id2 = Curves(i,2);
    New_v1 = NewPoints(id1,:);
    New_v2 = NewPoints(id2,:);
    for j = 1:nrcurves
        jd1 = Curves(j,1);
        jd2 = Curves(j,2);
        v1 = Points(jd1,:);
        v2 = Points(jd2,:);
        if (v1 == New_v1) 
            if (v2 == New_v2)
                flag = 1;
                break
            end
        elseif (v2 == New_v1) 
            if (v1 == New_v2)
                flag = 1;
                break
            end
        end
    end
    if (flag == 0)
        k = k+1;
        NewCurves(k,:) = Curves(i,:);
    end
end

NewCurves = NewCurves + nrpoints;

Curves = [Curves ; NewCurves];
Points = [Points ; NewPoints];

%Modify indices of new curves according to point to the
%original points in case of duplicates
for i = nrcurves+1:nrcurves+k
    id1 = Curves(i,1);
    id2 = Curves(i,2);
    New_v1 = Points(id1,:);
    New_v2 = Points(id2,:);
    for j = 1:nrpoints
        v = Points(j,:);
        if (v == New_v1)
            Curves(i,1) = j;
         elseif (v == New_v2)
            Curves(i,2) = j;
         end
    end
end
% % 

figure(2)
% nrcurves = nrcurves+k;
for i = 1:nrcurves
    id1 = Curves(i,1);
    id2 = Curves(i,2);
    v1 = Points(id1,:);
    v2 = Points(id2,:);
    v =  [v1;v2];
    plot3(v(:,1), v(:,2), v(:,3), '-sk','LineWidth',2);
    hold on
end
for i = nrcurves+1:nrcurves+k
    id1 = Curves(i,1);
    id2 = Curves(i,2);
    v1 = Points(id1,:);
    v2 = Points(id2,:);
    v =  [v1;v2];
    plot3(v(:,1), v(:,2), v(:,3), '-or','LineWidth',2);
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
% % 
% % 


