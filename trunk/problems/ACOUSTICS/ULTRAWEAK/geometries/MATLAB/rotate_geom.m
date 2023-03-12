function [nrpoints, Points, ...
          nrcurves,nrlines,nrcurvlines, Curves, ...
          nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
          nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume] ...
          = rotate_geom(nrpoints, Points, nrcurves,nrlines,nrcurvlines, Curves, ...
            nrtria,nrquad,nrrect,nrtraqua,nrptirec, Surfaces, ...
            nrprism,nrhexas,nrlinhex,nrtrahex, nrtetra, nrpyram, Volume)

% % Rotation    
% theta(1)= 0; 
% theta(2)= 0;
% theta(3)=90;
% rad   = theta*pi/180;
% sina(1:3) = sin(rad(1:3));
% cosa(1:3) = cos(rad(1:3));
%     
% % Rz = [cosa(3) -sina(3) 0 ; sina(3) cosa(3) 0 ; 0 0 1];
Rz = [0 -1 0 ; 1 0 0 ; 0 0 1];  

% Rz = [1 0 0 ; 0 1 0 ; 0 0 1];  

nrpoints_new = nrpoints;
for i = 1: nrpoints
    newpoint = (Rz*Points(i,:)')';
    %Check if the point already exists
    [Lia,Locb] = ismember(newpoint,Points, 'rows');
    if Lia
       inj(i) = Locb;
    else
       nrpoints_new = nrpoints_new+1;
       Points = [Points ; newpoint];
       inj(i) = nrpoints_new;
    end
end

%Construct new curves
nrcurves_new= nrcurves;
nrlines_new = nrlines;
nrcurvlines_new = nrcurvlines;

Curves_lines = Curves(1:nrlines,:);
Curves_curvlines = Curves(nrlines+1:end,:);
for i =1:nrlines
    id1 = Curves_lines(i,1);
    id2 = Curves_lines(i,2);
    newid1 = inj(id1);
    newid2 = inj(id2);
    jd1  = [newid1 newid2];
    flag = 0;
    for j = 1:nrlines_new
        int = intersect(jd1, Curves_lines(j,:));
        if (length(int) == 2)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
%total number of curves
       nrcurves_new = nrcurves_new + 1;
       nrlines_new = nrlines_new + 1;
%total number of lines and curvelines       
       newid = jd1;
       Curves_lines = [Curves_lines ; newid];
    end
end    

for i =1:nrcurvlines
    id1 = Curves_curvlines(i,1);
    id2 = Curves_curvlines(i,2);
    newid1 = inj(id1);
    newid2 = inj(id2);
    jd1  = [newid1 newid2];
    flag = 0;
    for j = 1:nrcurvlines_new
        int = intersect(jd1, Curves_curvlines(j,:));
        if (length(int) == 2)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
%total number of curves
       nrcurves_new = nrcurves_new + 1;
       nrcurvlines_new = nrcurvlines_new + 1;
%total number of lines and curvelines       
       newid = jd1;
       Curves_curvlines = [Curves_curvlines ; newid];
    end
end  


Curves = [Curves_lines; Curves_curvlines];


nrpoints = nrpoints_new;
nrlines =  nrlines_new;
nrcurvlines = nrcurvlines_new;
nrcurves = nrcurves_new;




















%Construct new surfaces
nrquad_new= nrquad;
nrrect_new = nrrect;
nrtraqua_new = nrtraqua;
nrptirec_new = nrptirec;

Surfaces_rect = Surfaces(1:nrrect,:);
Surfaces_traqua = Surfaces(nrrect+1:nrrect+nrtraqua,:);
Surfaces_ptirec = Surfaces(nrrect+nrtraqua+1:nrquad,:);

for i =1:nrrect
    id1 = Surfaces_rect(i,1);
    id2 = Surfaces_rect(i,2);
    id3 = Surfaces_rect(i,3);
    id4 = Surfaces_rect(i,4);
    newid1 = inj(id1);
    newid2 = inj(id2);
    newid3 = inj(id3);
    newid4 = inj(id4);
    
    jd1  = [newid1 newid2 newid3 newid4];
    flag = 0;
    for j = 1:nrrect_new
        int = intersect(jd1, Surfaces_rect(j,:));
        if (length(int) == 4)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
%total number of curves
       nrquad_new = nrquad_new + 1;
       nrrect_new = nrrect_new+1;
%total number of lines and curvelines       
       newid = jd1;
       Surfaces_rect = [Surfaces_rect ; newid];
    end
end
for i =1:nrtraqua
    id1 = Surfaces_traqua(i,1);
    id2 = Surfaces_traqua(i,2);
    id3 = Surfaces_traqua(i,3);
    id4 = Surfaces_traqua(i,4);
    newid1 = inj(id1);
    newid2 = inj(id2);
    newid3 = inj(id3);
    newid4 = inj(id4);
    
    jd1  = [newid1 newid2 newid3 newid4];
    flag = 0;
    for j = 1:nrtraqua_new
        int = intersect(jd1, Surfaces_traqua(j,:));
        if (length(int) == 4)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
%total number of curves
       nrquad_new = nrquad_new + 1;
       nrtraqua_new = nrtraqua_new+1;
%total number of lines and curvelines       
       newid = jd1;
       Surfaces_traqua = [Surfaces_traqua ; newid];
    end
end

for i =1:nrptirec
    id1 = Surfaces_ptirec(i,1);
    id2 = Surfaces_ptirec(i,2);
    id3 = Surfaces_ptirec(i,3);
    id4 = Surfaces_ptirec(i,4);
    newid1 = inj(id1);
    newid2 = inj(id2);
    newid3 = inj(id3);
    newid4 = inj(id4);
    
    jd1  = [newid1 newid2 newid3 newid4];
    flag = 0;
    for j = 1:nrptirec_new
        int = intersect(jd1, Surfaces_ptirec(j,:));
        if (length(int) == 4)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
%total number of curves
       nrquad_new = nrquad_new + 1;
       nrptirec_new = nrptirec_new+1;
%total number of lines and curvelines       
       newid = jd1;
       Surfaces_ptirec = [Surfaces_ptirec ; newid];
    end
end
Surfaces = [Surfaces_rect ; Surfaces_traqua ; Surfaces_ptirec];

nrquad= nrquad_new;
nrrect = nrrect_new;
nrtraqua = nrtraqua_new;
nrptirec = nrptirec_new;


%Construct new Polyhedra
 nrprism_new  = nrprism;
 nrhexas_new  = nrhexas;
 nrlinhex_new = nrlinhex;
 nrtrahex_new = nrtrahex;
 nrtetra_new  = nrtetra;
 nrpyram_new  = nrpyram;

Volume_linhex = Volume(1:nrlinhex, :);
Volume_trahex = Volume(nrlinhex+1:end, :);
for i =1:nrlinhex
    id1 = Volume_linhex(i,1);
    id2 = Volume_linhex(i,2);
    id3 = Volume_linhex(i,3);
    id4 = Volume_linhex(i,4);
    id5 = Volume_linhex(i,5);
    id6 = Volume_linhex(i,6);
    id7 = Volume_linhex(i,7);
    id8 = Volume_linhex(i,8);
    newid1 = inj(id1);
    newid2 = inj(id2);
    newid3 = inj(id3);
    newid4 = inj(id4);
    newid5 = inj(id5);
    newid6 = inj(id6);
    newid7 = inj(id7);
    newid8 = inj(id8);
    
    jd1  = [newid1 newid2 newid3 newid4 newid5 newid6 newid7 newid8];
    flag = 0;
    for j = 1:nrlinhex_new
        int = intersect(jd1, Volume_linhex(j,:));
        if (length(int) == 8)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
       nrhexas_new = nrhexas_new + 1;
       nrlinhex_new = nrlinhex_new + 1;
       newid = jd1;
       Volume_linhex  = [Volume_linhex ; newid];
    end
end    
for i =1:nrtrahex
    id1 = Volume_trahex(i,1);
    id2 = Volume_trahex(i,2);
    id3 = Volume_trahex(i,3);
    id4 = Volume_trahex(i,4);
    id5 = Volume_trahex(i,5);
    id6 = Volume_trahex(i,6);
    id7 = Volume_trahex(i,7);
    id8 = Volume_trahex(i,8);
    newid1 = inj(id1);
    newid2 = inj(id2);
    newid3 = inj(id3);
    newid4 = inj(id4);
    newid5 = inj(id5);
    newid6 = inj(id6);
    newid7 = inj(id7);
    newid8 = inj(id8);
    
    jd1  = [newid1 newid2 newid3 newid4 newid5 newid6 newid7 newid8];
    flag = 0;
    for j = 1:nrtrahex_new
        int = intersect(jd1, Volume_trahex(j,:));
        if (length(int) == 8)
            flag = 1;
            break
        end
    end    
     if (flag == 0)
       nrhexas_new = nrhexas_new + 1;
       nrtrahex_new = nrtrahex_new + 1;
       newid = jd1;
       Volume_trahex  = [Volume_trahex ; newid];
    end
end    


Volume = [Volume_linhex ; Volume_trahex];

% [nrhexas, nrlinhex, nrtrahex]


 nrprism  = nrprism_new;
 nrhexas  = nrhexas_new;
 nrlinhex = nrlinhex_new;
 nrtrahex = nrtrahex_new;
 nrtetra  = nrtetra_new;
 nrpyram  = nrpyram_new;

%  [nrhexas, nrlinhex, nrtrahex]

