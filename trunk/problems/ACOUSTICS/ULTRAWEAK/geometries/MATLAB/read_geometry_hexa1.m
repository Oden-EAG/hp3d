clear ; clc ;close all
format long

% read file
fid = fopen('hexa1');
FC = textscan(fid, '%s');
a = FC{1};
fclose(fid);


% 
% [nrpoints,nrcurves,nrlines,nrcurvlines,   ...
%  nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
%  nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, ...
%  Points,Curves,Surfaces,Volume] = read_file(a);
%     

nrpoints = a(9);
nrpoints = str2num(nrpoints{1});
b = a(11:end);
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
nrcurvlines = nrcurves-nrlines;


%------------------------------------------------
%               READ SURFACES
%------------------------------------------------
% Triangles
% j = j+4;
j=1;
nrtria = d(j); 
nrtria = str2num(nrtria{1});
j = j+2;
e = d(j:end);

nrquad = e(1);
nrquad = str2num(nrquad{1});
j = 5;
%Rectangles first
for i = 1:nrquad
    for k = 1:4
        y = e(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Surfaces(i,:) = x;
    str = e(j+5);
    str = str{1};
    if isequal(str, 'BilQua')
        j = j+7;
    else
        j = j+7;
      break 
    end    
end
nrrect = i;
%
%TraQuad 
ii = 0;
for i = nrrect+1:nrquad
    ii=ii+1;
    for k = 1:4
        y = e(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Surfaces(i,:) = x;
    str = e(j+5);
    str = str{1};
    if isequal(str, 'TraQua')
        j = j+7;
    else
        j = j+7;
      break 
    end    
end
nrtraqua = ii;
%TraQuad 
ii=0;
for i = nrrect+nrtraqua+1:nrquad
    ii=ii+1;
    for k = 1:4
        y = e(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Surfaces(i,:) = x;
    str = e(j+6);
    str = str{1};
    if isequal(str, 'PTIRec')
        j = j+8;
    else
        j = j+6;
      break 
    end    
end
nrptirec = ii;
j=j-2;
%------------------------------------------------
%               READ VOLUME SHAPES
%------------------------------------------------
%PRISMS
nrprism = e(j);
nrprism = str2num(nrprism{1});
j = j+2;

%HEXAS
g = e(j:end);
nrhexas = g(1);
nrhexas = str2num(nrhexas{1});

j = 6;
for i = 1:nrhexas
    for k = 1:8
        y = g(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Volume(i,:) = x;
    str = g(j+9);
    str = str{1};
    if isequal(str, 'Linear')
        j = j+12;
    else
        j = j+12;
      break 
    end    
end
nrlinhex = i;

%TraHex
ii=0;
for i = nrlinhex+1:nrhexas
    ii=ii+1;
    for k = 1:8
        y = g(j+k);
        y = y{1};
        x(k) = str2num(y);
    end
    Volume(i,:) = x;
    str = g(j+9);
    str = str{1};
    if isequal(str, 'TraHex')
        j = j+12;
    else
        j = j+9;
      break 
    end    
end
nrtrahex = ii;
j=j-3;
nrtetra = g(j);
nrtetra = str2num(nrtetra{1});
j = j+2;

nrpyram = g(j);
nrpyram = str2num(nrpyram{1});


disp('translate in x')
for i = 1:6
    i
[nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
 = translate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
               nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
               nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume,1,0,0);

end

disp('translate in y')
for i = 1:6
    i
[nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
 = translate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
               nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
               nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume,0,1,0);

end
% 

disp('translate in z')
for i = 1:6
    i
[nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
 = translate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
               nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
               nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume,0,0,1);

end
% 
disp('SUMMARY')

%Points
fprintf('Number of Points     = %i \n', nrpoints);

%Curves
fprintf('Number of Curves     = %i \n', nrcurves);
fprintf('Number of Seglin     = %i \n', nrlines);
fprintf('Number of SurfCur    = %i \n', nrcurvlines);

%Surfaces]
fprintf('Number of Triangles  = %i \n', nrtria);
fprintf('Number of Quads      = %i \n', nrquad);
fprintf('Number of BilQua     = %i \n', nrrect);
fprintf('Number of TraQua     = %i \n', nrtraqua);
fprintf('Number of PTIRec     = %i \n', nrptirec);


%Volume
fprintf('Number of Prisms     = %i \n', nrprism);
fprintf('Number of Hexas      = %i \n', nrhexas);
fprintf('Number of Linear hex = %i \n', nrlinhex);
fprintf('Number of Trahex     = %i \n', nrtrahex);



%remove 365rd hexa (the middle one)
nrlinhex = nrlinhex-1;
nrhexas = nrhexas-1;
Volume = [Volume(1:171,:) ; Volume(173:end,:)];
Points = Points/7;
plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
filename = 'hexa342_orient0_unit_cube';
write_file(filename,nrpoints,nrcurves,nrlines,nrcurvlines,   ...
           nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
           nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
           Points,Curves,Surfaces,Volume);


