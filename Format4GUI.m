function [nt1,nt2,nt3]=Format4GUI(neuront1,neuront2,neuront3,reg,keeps,vfall,trinum1,trinum2,trinum3);

t1k=keeps{1};
t2k=keeps{2};
t3k=keeps{3};

vf1=vfall{1};
vf2=vfall{2};
vf3=vfall{3};

neuront1k.c=neuront1.C(:,t1k)'+neuront1.YrA(:,t1k)';
neuront2k.c=neuront2.C(:,t2k)'+neuront2.YrA(:,t2k)';

method='notrial';
neuronst1=GetCS(neuront1k.c,vf1); 
neuronst2=GetCS(neuront2k.c,vf2); 

nt1.c=neuronst1{1};
nt1.s=neuronst1{2};
nt1.center=neuront1.coor(t1k,:);
nt1.center=[nt1.center(:,2) nt1.center(:,1)];
nt1.Cn=reshape(neuront1.all_A,neuront1.dims(1),neuront1.dims(2));

nt2.c=neuronst2{1};
nt2.s=neuronst2{2};
nt2.center=neuront2.coor(t2k,:);
nt2.center=[nt2.center(:,2) nt2.center(:,1)];
nt2.Cn=reshape(neuront2.all_A,neuront2.dims(1),neuront2.dims(2));

neuront3k.c=neuront3.C(:,t3k)'+neuront3.YrA(:,t3k)';
neuronst3=GetCS(neuront3k.c,vf3); 
nt3.c=neuronst3{1};
nt3.s=neuronst3{2};
nt3.center=neuront3.coor(t3k,:);
nt3.center=[nt3.center(:,2) nt3.center(:,1)];
nt3.Cn=reshape(neuront3.all_A,neuront3.dims(1),neuront3.dims(2));

end

