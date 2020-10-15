% close all
clear all

% name
prob='quarter_sphere_6';

sfx='';
% does this problem has interior volume seeds? true or false
intvol=true;

inc_factor=1;

% Determine location of VoroCrust ouput files
folder='/org/groups/aeg/jaime_mora/vorocrust_repo/vorocrust/examples/';
folder=[folder,prob,'/output'];
% read data
intvol_n=0;
if intvol
    seeds_intvol=csvread(strcat(folder,'/Interior_volume_seeds',sfx,'.csv'),1,0);
    intvol_n=size(seeds_intvol,1);
end
seeds_int=csvread(strcat(folder,'/Interior_seeds',sfx,'.csv'),1,0);
seeds_ext=csvread(strcat(folder,'/Exterior_seeds',sfx,'.csv'),1,0);
seeds_extvol=csvread(strcat(folder,'/Exterior_volume_seeds',sfx,'.csv'),1,0);
% save number of cells of each type
% 
int_n=size(seeds_int,1);
ext_n=size(seeds_ext,1);
extvol_n=size(seeds_extvol,1);


if intvol
    vorocrust_seeds=[seeds_intvol(:,:);
                     seeds_int(:,:);
                     seeds_ext(:,:);
                     seeds_extvol(:,:)];
else
    vorocrust_seeds=[seeds_int(:,:);
                     seeds_ext(:,:);
                     seeds_extvol(:,:)];
end
cells_n=intvol_n+int_n;
h_K=zeros(cells_n,1);ratio_K=h_K;
[V_vc,C_vc]=voronoin(vorocrust_seeds(1:cells_n+ext_n,1:3));
for k=1:cells_n
    poly=C_vc{k};
    poly_v=V_vc(poly,:);
    dmax=0;
    for i=1:length(poly)
        for j=1:length(poly)
            dij=norm(poly_v(i,:)-poly_v(j,:));
            dmax=max(dmax,dij);
        end
    end
    h_K(k)=dmax;
    ratio_K(k)=(dmax/vorocrust_seeds(k,4));
end
clear V_vc C_vc





elem_ref_flag=readmatrix(strcat('./output/elem_ref_flag',sfx));
total_elem = size(elem_ref_flag,1);
elem_ref_factor=ones(total_elem,1)-0.5*elem_ref_flag(:,2);
% 
% el=0;
% 
% if intvol
%     new_intvol=zeros(intvol_n,4);
%     rmin_intvol=min(seeds_intvol(:,4));
%     rmax_intvol=max(seeds_intvol(:,4));
%     r_threshold=rmin_intvol+0.75*(rmax_intvol-rmin_intvol);
%     for j=1:intvol_n
%         el = el + 1;
%         new_intvol(j,1:3)=seeds_intvol(j,1:3);
% %         new_intvol(j,4)=seeds_intvol(j,4)*elem_ref_factor(el);
%     if seeds_intvol(j,4) >= r_threshold
%         if elem_ref_flag(el,2)==1
% %             new_int(j,4)=seeds_int(j,4)*0.5*inc_factor;
%             new_intvol(j,4)=rmax_intvol*0.5*inc_factor;
%         else
%             new_intvol(j,4)=rmax_intvol*inc_factor;
%         end         
% %         new_intvol(j,4)=seeds_intvol(j,4)*elem_ref_factor(el);
%     else
%         if elem_ref_flag(el,2)==1
% %             new_intvol(j,4)=seeds_intvol(j,4)*0.5*inc_factor;
%             new_intvol(j,4)=0.5*(seeds_intvol(j,4)+rmax_intvol)*0.5*inc_factor;
%         else
%             new_intvol(j,4)=0.5*(seeds_intvol(j,4)+rmax_intvol)*inc_factor;
%         end 
%     end          
%     end
% end
% 
% 
% new_int=zeros(int_n,4);
% % get the minimum radius
% rmin_surf=min(seeds_int(:,4));
% rmax_surf=max(seeds_int(:,4));
% r_threshold=rmin_surf+0.75*(rmax_surf-rmin_surf);
% 
% for j=1:int_n
%     el = el + 1;
%     new_int(j,1:3)=seeds_int(j,1:3);
%     % Usually surface seeds are much smaller than edge and corner seeds.
%     % To prevent from having too small spheres for the surfaces, 
%     % if the current sphere is smaller than the minimum multiplied by a
%     % factor between 1 and 2, we do not refine. 
%     if seeds_int(j,4) >= r_threshold
%         if elem_ref_flag(el,2)==1
% %             new_int(j,4)=seeds_int(j,4)*0.5*inc_factor;
%             new_int(j,4)=rmax_surf*0.5*inc_factor;
%         else
%             new_int(j,4)=rmax_surf*inc_factor;
%         end         
% %         new_int(j,4)=seeds_int(j,4)*elem_ref_factor(el);
%     else
%         if elem_ref_flag(el,2)==1
% %             new_int(j,4)=seeds_int(j,4)*0.5*inc_factor;
%             new_int(j,4)=0.5*(seeds_int(j,4)+rmax_surf)*0.5*inc_factor;
%         else
%             new_int(j,4)=0.5*(seeds_int(j,4)+rmax_surf)*inc_factor;
%         end 
%     end
% end

el=0;

if intvol
    h_K_intvol=h_K(1:intvol_n);
    new_intvol=zeros(intvol_n,4);
    hmin_intvol=min(h_K_intvol);
    hmax_intvol=max(h_K_intvol);
    h_threshold=hmin_intvol+0.5*(hmax_intvol-hmin_intvol);
    for j=1:intvol_n
        el = el + 1;
        new_intvol(j,1:3)=seeds_intvol(j,1:3);
%         new_intvol(j,4)=h_K_intvol(j)*max(0.5,1-elem_ref_flag(el,3)^0.5);
        if h_K_intvol(j)>=h_threshold
            if elem_ref_flag(el,2)==1
                new_intvol(j,4)=(h_K_intvol(j))*0.5*inc_factor;
            else
                new_intvol(j,4)=(h_K_intvol(j))*inc_factor;
            end         
        else
            if elem_ref_flag(el,2)==1
                new_intvol(j,4)=(h_K_intvol(j)+h_threshold)*0.5*inc_factor;
            else
                new_intvol(j,4)=0.5*(h_K_intvol(j)+h_threshold)*inc_factor;
            end         
        end
%     
    end
end

h_K_int=h_K(1:int_n);
new_int=zeros(int_n,4);
% get the minimum radius
hmin_surf=min(h_K_int);
hmax_surf=max(h_K_int);
h_threshold=hmin_surf+0.5*(hmax_surf-hmin_surf);

for j=1:int_n
    el = el + 1;
    new_int(j,1:3)=seeds_int(j,1:3);
%     new_int(j,4)=h_K_int(j)*max(0.5,1-elem_ref_flag(el,3)^0.5);
    % Usually surface seeds are much smaller than edge and corner seeds.
    % To prevent from having too small spheres for the surfaces, 
    % if the current sphere is smaller than the minimum multiplied by a
    % factor between 1 and 2, we do not refine. 
        if h_K_int(j)>=h_threshold
            if elem_ref_flag(el,2)==1
                new_int(j,4)=(h_K_int(j))*0.5*inc_factor;
            else
                new_int(j,4)=(h_K_int(j))*inc_factor;
            end         
        else
            if elem_ref_flag(el,2)==1
                new_int(j,4)=0.5*(h_K_int(j)+h_threshold)*0.5*inc_factor;
            else
                new_int(j,4)=0.5*(h_K_int(j)+h_threshold)*inc_factor;
            end         
        end        
end

new_ext(:,1:3)=seeds_ext(:,1:3);
new_ext(:,4)=new_int(:,4);

% 
% if intvol
%     file0_label=strcat(folder,'/Interior_volume_NEW_SIZE.csv');
%     file0=fopen(file0_label,'w');
%     for jv=1:intvol_n
%         fprintf(file0,'%17.15f,%17.15f,%17.15f,%17.15f\n',new_intvol(jv,1:4));
%     end
% end
% 
% 
% file1_label=strcat(folder,'/Interior_NEW_SIZE.csv');
% file1=fopen(file1_label,'w');
% for jv=1:int_n
%    fprintf(file1,'%17.15f,%17.15f,%17.15f,%17.15f\n',new_int(jv,1:4));
% end

if intvol
    all_spheres=[new_intvol;new_int;new_ext];
else
    all_spheres=[new_int;new_ext];
end
n_all=size(all_spheres,1);
file1_label=strcat(folder,'/NEW_SIZE.csv');
file1=fopen(file1_label,'w');
for jv=1:n_all
   fprintf(file1,'%17.15f\t%17.15f\t%17.15f\t%17.15f\n',all_spheres(jv,1:4));
end




[V_vc,C_vc]=voronoin(all_spheres(:,1:3));
% construct all facets in interior volume and interior surface
V_vor=[];
K=[];
nvf=[];    
nfacets=0;
nvtotal=0;
lim_nvf=20;
[V_vc,C_vc]=voronoin(all_spheres(:,1:3));


first_cell=1;
if intvol
    last_cell=intvol_n+int_n
%     rmin=min(rmin_intvol,rmin_surf);
%     rmax=max(rmax_intvol,rmax_surf);
    rmin=min(hmin_intvol,hmin_surf);
    rmax=max(hmax_intvol,hmax_surf);
else
    last_cell=int_n
%     rmin=rmin_surf;
%     rmax=rmax_surf;    
    rmin=hmin_surf;
    rmax=hmax_surf; 
end


% free some memory
clear seeds_intvol seeds_int seeds_ext NFC CELLS F

ncells_plot=last_cell-(first_cell-1);
for i=first_cell:last_cell
    verts_in_cell=[];
    verts_in_cell=V_vc(C_vc{i},:);
%     if sum(find(~(C_vc{i}-812)))~=0 
%         disp(i)
%     end
%     if i==539
%         verts_in_cell
%         C_vc{i}
%     end
    nvc=size(verts_in_cell,1);
    V_vor=[V_vor;verts_in_cell]; 
%     nv=size(V_vc(C_vc{i},:),1);
    Kmin=minConvexHull(verts_in_cell);
    nk=length(Kmin);
    Ktmp=zeros(nk,lim_nvf);
    for j=1:nk
        nvf(j+nfacets)=length(Kmin{j});
        Ktmp(j,1:nvf(j+nfacets))=Kmin{j}+nvtotal;
    end
    K=[K; Ktmp];
    CELLS{i}=nfacets+1:nfacets+nk;
    NFC(i)=nk;
    nfacets=nfacets+nk;
    nvtotal = nvtotal + nvc;
end
max_nvf=max(nvf);
F=K(:,1:max_nvf);
clear K;

[V_uni,ia_v,ic_v]=unique(V_vor,'rows','stable');

% F_uni=NaN*ones(nfu,max_nvf);
for i=1:nfacets
    F(i,1:nvf(i))=ic_v(F(i,1:nvf(i)));
%     for j=1:max_nvf
%         if F_uni(i,j)==Inf
%             F_uni(i,j)=NaN;
%         end
%     end
end   


[~,ia_f,ic_f]=unique(sort(F,2),'rows','stable');
nfu=length(ia_f);

F_uni=NaN*ones(nfu,max_nvf);
for i=1:nfu
    nvf_uni(i)=nvf(ia_f(i));
    F_uni(i,1:nvf_uni(i))=F(ia_f(i),1:nvf_uni(i))';
end


FACE_BC=2*ones(nfu,1);
for jc=1:ncells_plot
    facets=[];
    facets=CELLS{jc};
    CELLS{jc}=ic_f(facets)';
    CELL_F_ORI{jc}=zeros(1,length(facets));
    for jf=1:length(facets)
%         disp(jf)
        kf=ic_f(facets(jf));
        FACE_BC(kf)=FACE_BC(kf)-1;
        if FACE_BC(kf)==0
            CELL_F_ORI{jc}(jf)=1;
        end
    end
end

% colorf=zeros(size(F,1),3);
% for jc=1:ncells_plot
%     facets=[];
%     facets=CELLS{jc};
%     cellcolor=rand(1,3);
%     for jf=1:NFC(jc)
%         colorf(facets(jf),1:3)=cellcolor;
%     end
% end
% 
% patch('Faces',F,'Vertices',V_vor,'FaceVertexCData',colorf,'FaceColor','flat','FaceAlpha',0.85)

colorf=zeros(size(F_uni,1),3);
for jc=1:ncells_plot
    facets=[];
    facets=CELLS{jc};
    if elem_ref_flag(jc,2)==1
        cellcolor=[1 0 0];%rand(1,3)*0.3+.7;
    else
        cellcolor=[0.95 0.95 0.95];
    end
%     cellcolor=[0.95 0.95*(1-elem_ref_flag(jc,3)^0.5) 0.95*(1-elem_ref_flag(jc,3)^0.5)];
%     cellcolor2=[0.95 0.95*(all_spheres(jc,4)/rmax) 0.95*(all_spheres(jc,4)/rmax)];
    for jf=1:NFC(jc)
        colorf(facets(jf),1:3)=cellcolor;
%         colorf2(facets(jf),1:3)=cellcolor2;
    end
end
% TEMP=CELLS{539};
figure(8)
patch('Faces',F_uni(:,:),'Vertices',V_uni,'FaceVertexCData',colorf(:,:),'FaceColor','flat','FaceAlpha',0.8)
grid on
daspect([1 1 1])
view(23,23)
lightangle(20,60)
% figure(15)
% patch('Faces',F_uni(:,:),'Vertices',V_uni,'FaceVertexCData',colorf2(:,:),'FaceColor','flat','FaceAlpha',0.8)
% grid on
% daspect([1 1 1])
% view(23,23)
% lightangle(20,60)