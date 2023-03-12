function write_file(fname,nrpoints,nrcurves,nrlines,nrcurvlines,   ...
 nrtria,nrquad,nrrect,nrtraqua,nrptirec, ...
 nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram,  ...
 Points,Curves,Surfaces,Volume)


% %check if the curve numbering is correct
% if (nrlines+nrcurvlines /= nrcurves) 
%     disp('Inconsistency with number of curves')
%     fprintf('nrlines = %i, nrcurvlines = %i, nrcurves = %i', nrlines, nrcurvlines, nrcurves)
% end
% %check if the surface numbering is correct
% if (nrrect+nrtraquad+nrtrahex /= nrquad) 
%     disp('Inconsistency with number of surfaces')
%     fprintf('nrrect = %i, nrtraquad = %i, nrtrahex = %i, nrquad = %i', nrrect, nrtraquad, nrtrahex, nrquad)
% end


filename = fopen(fname,'w');
fprintf(filename,'3  3 NDIM, MANDIM \n\n');
fprintf(filename,'1    no surfaces \n\n');
fprintf(filename,['Sphere    surface 1: sphere \n', ...
                 '0.d0 0.d0 2.0d0     sphere center \n', ...
                 '0.25d0              sphere radius \n\n']);
fprintf(filename,'1         NRDOMAIN \n\n');

%POINTS
fprintf(filename,'%i        NRPOINT  \n\n', nrpoints);
for i = 1: nrpoints
    fprintf(filename, 'Regular             point %i \n ', i);
    fprintf(filename, '%.16fd0 %.16fd0 %.16fd0  \n\n' , Points(i,1:3));
end



%CURVES
fprintf(filename,'%i        NRCURVE  \n\n', nrcurves);

for i = 1: nrlines
    fprintf(filename, 'Seglin             curve %i \n ', i);
    fprintf(filename, '%i %i \n\n' , Curves(i,1:2));
end
for i = nrlines+1: nrcurves
    fprintf(filename, '1SurfsCur             curve %i \n ', i);
    fprintf(filename, '%i %i \n 1 \n\n' , Curves(i,1:2));
    
end

%SURFACES

fprintf(filename,'%i         NRTRIAN  \n\n', nrtria);
fprintf(filename,'%i        NRRECTA  \n\n', nrquad);

for i = 1: nrrect
    fprintf(filename, 'BilQua            rectangle %i \n ', i);
    fprintf(filename, '%i %i %i %i \n\n' , Surfaces(i,1:4));
end
for i = nrrect+1:nrrect+nrtraqua
    fprintf(filename, 'TraQua            rectangle %i \n ', i);
    fprintf(filename, '%i %i %i %i \n\n' , Surfaces(i,1:4));
end
for i = nrrect+nrtraqua+1 : nrquad
    fprintf(filename, 'PTIRec            rectangle %i \n ', i);
    fprintf(filename, '%i %i %i %i \n 1 \n\n' , Surfaces(i,1:4));
end

%VOLUME
fprintf(filename,'%i         NRPRISM  \n\n', nrprism);
fprintf(filename,'%i        NRHEXAS  \n\n', nrhexas);

for i = 1: nrlinhex
    fprintf(filename, 'Linear            hexa %i \n ', i);
    fprintf(filename, '1   %i %i %i %i %i %i %i %i \n\n' , Volume(i,1:8));
end
for i = nrlinhex+1 : nrhexas
    fprintf(filename, 'TraHex            hexa %i \n ', i);
    fprintf(filename, '1   %i %i %i %i %i %i  %i %i \n\n' , Volume(i,1:8));
end


fprintf(filename,'%i    NRTETRA  \n\n', nrtetra);
fprintf(filename,'%i    NRPYRAM  \n\n', nrpyram);

fprintf(filename, ['0    number of curves to rotate \n\n', ...  
                   '0    number of triangles to rotate \n\n', ...
                   '0    number of rectangles to rotate \n']);


fclose(filename);
% type hemisphere_matlab
