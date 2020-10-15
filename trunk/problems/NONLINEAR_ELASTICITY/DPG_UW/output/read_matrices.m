close all

gram=[];dg=[];
gram=csvread('gram');
gram(:,end)=[];
dg=diag(diag(gram));
gram=gram+gram'-dg;

gram_L=[];
gram_L=csvread('gram_fact');
gram_L(:,end)=[];

EnrEverything=[];
EnrEverything=csvread('EnrEverything');
EnrEverything(:,end)=[];

EnrEverything_post=[];
EnrEverything_post=csvread('EnrEverything_post');
EnrEverything_post(:,end)=[];

EnrStiffness = [];
EnrStiffness = EnrEverything(:,1:end-1);

FullDPG=[];
FullDPG=csvread('FullDPG');
FullDPG(:,end)=[];

FOSLS=[];
FOSLS=csvread('FOSLS');
FOSLS(:,end)=[];

gram_rec=gram_L*gram_L';
gramdiff=abs(gram-gram_rec);

gram_filt=zeros(189);
for i=1:189
    rows=[];
    rows=find(gramdiff(:,i)>1E-4);
    gram_filt(rows,i)=gramdiff(rows,i);
end
spy(gram_filt)

figure(2)
GIB=gram\EnrEverything;
K=EnrStiffness'*GIB;
dpg_diff=abs(K-FullDPG);
dpg_filt=zeros(54,55);
for i=1:55
    rows=[];
    rows=find(dpg_diff(:,i)>1E-4);
    dpg_filt(rows,i)=dpg_diff(rows,i);
end
spy(dpg_filt)