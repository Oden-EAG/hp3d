clear; clc; close all
format long
filename=strcat('hemisphere');
a=fileread(filename);


% read file
fid = fopen('hemisphere');
FC = textscan(fid, '%s');
a = FC{1};
fclose(fid);



[nrpoints,nrcurves,nrlines,nrcurvlines,   ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, ...
 Points,Curves,Surfaces,Volume] = read_file(a);
    


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

% plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
% write_file(nrpoints,nrcurves,nrlines,nrcurvlines,   ...
%            nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
%            nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
%            Points,Curves,Surfaces,Volume)

% 
% % plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
% filename = 'quartersphere_matlab';
% % write_file(filename,nrpoints,nrcurves,nrlines,nrcurvlines,   ...
% %            nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
% %            nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
% %            Points,Curves,Surfaces,Volume)
%        
% % plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
% 
% 
% 
% filename = 'hemisphere_matlab';
[nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
 = rotate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
               nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
               nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume);
% % write_file(filename,nrpoints,nrcurves,nrlines,nrcurvlines,   ...
% %            nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
% %            nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
% %            Points,Curves,Surfaces,Volume)
% 
% % plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
% 
% 
% 
% 
% 
% filename = '3_4sphere_matlab';
[nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
 = rotate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
               nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
               nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume);

% % write_file(filename,nrpoints,nrcurves,nrlines,nrcurvlines,   ...
% %            nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
% %            nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
% %            Points,Curves,Surfaces,Volume)
% 
% % plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
% 
% 
% 
filename = 'sphere_4unit_cube';
[nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
 = rotate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
               nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
               nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume);
Points(:,1:2) = Points(:,1:2) +2;
% Points = Points/4;
write_file(filename,nrpoints,nrcurves,nrlines,nrcurvlines,   ...
           nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
           nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
           Points,Curves,Surfaces,Volume)
% 
plot_geom(nrcurves, Curves,nrquad,Surfaces,nrhexas,Volume,Points);
% 
