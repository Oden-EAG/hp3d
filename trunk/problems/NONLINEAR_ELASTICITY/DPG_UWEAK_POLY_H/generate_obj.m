clear all
close all
n=24;
dt=pi/2/n;
phi=-pi/2:dt:pi/2;
theta=pi/2:-dt:0;
NUMV=16+(2*n-1)*(n+1);
VERTS=zeros(NUMV,3);
p=16;
a=15; b=16;

VERTS(1 ,:)=[0 0 -2];
VERTS(2 ,:)=[0 0 2];
VERTS(3 ,:)=[0 1 -1];
VERTS(4 ,:)=[0 1 1];
VERTS(5 ,:)=[0 2 -2];
VERTS(6 ,:)=[0 2 0];
VERTS(7 ,:)=[0 2 2];
VERTS(8 ,:)=[2 2 -2];
VERTS(9 ,:)=[2 2 2];
VERTS(10,:)=[2 0 -2];
VERTS(11,:)=[2 0 0];
VERTS(12,:)=[2 0 2];
VERTS(13,:)=[1 0 -1];
VERTS(14,:)=[1 0 1];
VERTS(15,:)=[0 0 -1];
VERTS(16,:)=[0 0 1];

NUMF=18+n+n+2*n*(2*n-2)+8+4*n;
FACES=zeros(NUMF,3);
f=18;


FACES(1 ,:)=[1 3 5];
FACES(2 ,:)=[3 6 5];
FACES(3 ,:)=[7 4 2];
FACES(4 ,:)=[6 4 7];
FACES(5 ,:)=[5 6 8];
FACES(6 ,:)=[8 6 9];
FACES(7 ,:)=[6 7 9];
FACES(8 ,:)=[2 9 7];
FACES(9 ,:)=[2 12 9];
FACES(10,:)=[1 5 8];
FACES(11,:)=[1 8 10];
FACES(12,:)=[8 9 11];
FACES(13,:)=[8 11 10];
FACES(14,:)=[9 12 11];
FACES(15,:)=[12 2 14];
FACES(16,:)=[11 12 14];
FACES(17,:)=[13 1 10];
FACES(18,:)=[11 13 10];



for j=1:n+1
    for i=2:2*n    
        p=p+1;
        [x, y, z]=sph2cart(theta(j),phi(i),1.0);
        VERTS(p,:)=[x y z];
        if (j==1) && i==n+1
            c=p;
        elseif (j==n+1) && i==n+1
            d=p;
        elseif (j==1) && i==n/2+1
            g=p;
        elseif (j==1) && i==3*n/2+1
            h=p;
        elseif (j==n+1) && i==n/2+1
            r=p;
        elseif (j==n+1) && i==3*n/2+1
            s=p;
        end
        if i==2 && j>1
            f=f+1;
            FACES(f,:)=[a p p-(2*n-1)];
        end
        if i==2*n && j>1
            f=f+1;
            FACES(f,:)=[b p-(2*n-1) p];
        end
        if j==1
           if i==2
               f=f+1;
               FACES(f,:)=[1 a p];
           elseif i<=n/2+1%p<=g
               f=f+1;
               FACES(f,:)=[p-1 p 1];
           elseif i>n/2+1 && i<=3*n/2+1%p>g && p<=h
               f=f+1;
               FACES(f,:)=[p-1 p 6];
           elseif i>3*n/2+1%p>h
               f=f+1;
               FACES(f,:)=[p-1 p 2];
               if i==2*n
                   f=f+1;
                   FACES(f,:)=[b 2 p];
               end
           end
        elseif j==n+1
            if i==2
               f=f+1;
               FACES(f,:)=[1 p a];
           elseif i<=n/2+1%p<=r
               f=f+1;
               FACES(f,:)=[p p-1 1];
           elseif i>n/2+1 && i<=3*n/2+1%p>r && p<=s
               f=f+1;
               FACES(f,:)=[p p-1 11];
           elseif i>3*n/2+1%p>s
               f=f+1;
               FACES(f,:)=[p p-1 2];
               if i==2*n
                  f=f+1;
                  FACES(f,:)=[b p 2];
               end
           end
        end
               
        if j>1 && i<2*n
            f=f+1;
            FACES(f,:)=[p-(2*n-1) p p-(2*n-1)+1];
            f=f+1;
            FACES(f,:)=[p-(2*n-1)+1 p p+1];
        end
    end
end

f=f+1;
FACES(f,:)=[1 g 3];
f=f+1;
FACES(f,:)=[3 g 6];
f=f+1;
FACES(f,:)=[6 h 4];
f=f+1;
FACES(f,:)=[4 h 2];
f=f+1;
FACES(f,:)=[1 13 r];
f=f+1;
FACES(f,:)=[r 13 11];
f=f+1;
FACES(f,:)=[11 14 s];
f=f+1;
FACES(f,:)=[s 14 2];
patch('Faces',FACES,'Vertices',VERTS,'FaceColor','Blue','FaceAlpha',0.85)
daspect([1 1 1])


filename=strcat('quarter_sphere_',num2str(n),'.obj');
file1=fopen(filename,'w');

fprintf(file1,strcat('#\t',filename,'\n'));
fprintf(file1,strcat('#\t',num2str(NUMV),' vertices\n'));
fprintf(file1,strcat('#\t',num2str(NUMF),' faces\n'));
fprintf(file1,'\n');

for jv=1:NUMV
    fprintf(file1,'v %17.15fd0\t%17.15fd0\t%17.15fd0\n',VERTS(jv,1:3));
end
for jf=1:NUMF
    fprintf(file1,'f %d\t%d\t%d\t\n',FACES(jf,1:3));
end