%% A comprehensive script for calcium imaging analysis of Social Memory tests
%Get CNMFE Data

filename1='CaIdCA2ID_L1N1_analysis_results.hdf5';
filename2='CaIdCA2ID_LL_analysis_results.hdf5';
filename3='CaIdCA2ID_NN_analysis_results.hdf5';

neuron_SRMt1=GetCNMFEData(filename1);
neuron_SRMt2=GetCNMFEData(filename2);
neuron_SRMt3=GetCNMFEData(filename3);


%% Combine Trials (3)

neuronall.C=[];
neuronall.S=[];
neuronall.YrA=[];

for n=1:length(cell_registered_struct.cell_to_index_map(:,1))
    if min(cell_registered_struct.cell_to_index_map(n,:))~=0
        neuronall.C(n,:)=[neuron_SRMt1.C(:,cell_registered_struct.cell_to_index_map(n,1)); neuron_SRMt2.C(:,cell_registered_struct.cell_to_index_map(n,2)); neuron_SRMt3.C(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.S(n,:)=[neuron_SRMt1.S(:,cell_registered_struct.cell_to_index_map(n,1)); neuron_SRMt2.S(:,cell_registered_struct.cell_to_index_map(n,2)); neuron_SRMt3.S(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.YrA(n,:)=[neuron_SRMt1.YrA(:,cell_registered_struct.cell_to_index_map(n,1)); neuron_SRMt2.YrA(:,cell_registered_struct.cell_to_index_map(n,2)); neuron_SRMt3.YrA(:,cell_registered_struct.cell_to_index_map(n,3))];
        
    elseif min(cell_registered_struct.cell_to_index_map(n,1:2))~=0
        neuronall.C(n,:)=[neuron_SRMt1.C(:,cell_registered_struct.cell_to_index_map(n,1)); neuron_SRMt2.C(:,cell_registered_struct.cell_to_index_map(n,2)); zeros(length(neuron_SRMt3.C(:,1)),1)];
        neuronall.S(n,:)=[neuron_SRMt1.S(:,cell_registered_struct.cell_to_index_map(n,1)); neuron_SRMt2.S(:,cell_registered_struct.cell_to_index_map(n,2)); zeros(length(neuron_SRMt3.C(:,1)),1)];
        neuronall.YrA(n,:)=[neuron_SRMt1.YrA(:,cell_registered_struct.cell_to_index_map(n,1)); neuron_SRMt2.YrA(:,cell_registered_struct.cell_to_index_map(n,2)); zeros(length(neuron_SRMt3.C(:,1)),1)];
    
    elseif min(cell_registered_struct.cell_to_index_map(n,2:3))~=0
        neuronall.C(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); neuron_SRMt2.C(:,cell_registered_struct.cell_to_index_map(n,2)); neuron_SRMt3.C(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.S(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); neuron_SRMt2.S(:,cell_registered_struct.cell_to_index_map(n,2)); neuron_SRMt3.S(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.YrA(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); neuron_SRMt2.YrA(:,cell_registered_struct.cell_to_index_map(n,2)); neuron_SRMt3.YrA(:,cell_registered_struct.cell_to_index_map(n,3))];

    elseif min(cell_registered_struct.cell_to_index_map(n,1))~=0&&min(cell_registered_struct.cell_to_index_map(n,3))~=0
        neuronall.C(n,:)=[neuron_SRMt1.C(:,cell_registered_struct.cell_to_index_map(n,1)); zeros(length(neuron_SRMt2.C(:,1)),1); neuron_SRMt3.C(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.S(n,:)=[neuron_SRMt1.S(:,cell_registered_struct.cell_to_index_map(n,1)); zeros(length(neuron_SRMt2.C(:,1)),1); neuron_SRMt3.S(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.YrA(n,:)=[neuron_SRMt1.YrA(:,cell_registered_struct.cell_to_index_map(n,1)); zeros(length(neuron_SRMt2.C(:,1)),1); neuron_SRMt3.YrA(:,cell_registered_struct.cell_to_index_map(n,3))];

    elseif max(cell_registered_struct.cell_to_index_map(n,1:2))==0
        neuronall.C(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); zeros(length(neuron_SRMt2.C(:,1)),1); neuron_SRMt3.C(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.S(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); zeros(length(neuron_SRMt2.C(:,1)),1); neuron_SRMt3.S(:,cell_registered_struct.cell_to_index_map(n,3))];
        neuronall.YrA(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); zeros(length(neuron_SRMt2.C(:,1)),1); neuron_SRMt3.YrA(:,cell_registered_struct.cell_to_index_map(n,3))];

    elseif max(cell_registered_struct.cell_to_index_map(n,2:3))==0
        neuronall.C(n,:)=[neuron_SRMt1.C(:,cell_registered_struct.cell_to_index_map(n,1)); zeros(length(neuron_SRMt2.C(:,1)),1); zeros(length(neuron_SRMt3.C(:,1)),1)];
        neuronall.S(n,:)=[neuron_SRMt1.S(:,cell_registered_struct.cell_to_index_map(n,1)); zeros(length(neuron_SRMt2.C(:,1)),1); zeros(length(neuron_SRMt3.C(:,1)),1)];
        neuronall.YrA(n,:)=[neuron_SRMt1.YrA(:,cell_registered_struct.cell_to_index_map(n,1)); zeros(length(neuron_SRMt2.C(:,1)),1); zeros(length(neuron_SRMt3.C(:,1)),1)];

    else
        neuronall.C(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); neuron_SRMt2.C(:,cell_registered_struct.cell_to_index_map(n,2)); zeros(length(neuron_SRMt3.C(:,1)),1)];
        neuronall.S(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); neuron_SRMt2.S(:,cell_registered_struct.cell_to_index_map(n,2)); zeros(length(neuron_SRMt3.C(:,1)),1)];
        neuronall.YrA(n,:)=[zeros(length(neuron_SRMt1.C(:,1)),1); neuron_SRMt2.YrA(:,cell_registered_struct.cell_to_index_map(n,2)); zeros(length(neuron_SRMt3.C(:,1)),1)];

    end
end

SRMt2keep=cell_registered_struct.cell_to_index_map(:,1);
SRMt2keepmod=SRMt2keep(SRMt2keep~=0);
SRMt2keepmodSRMt1=SRMt2keepmod(SRMt2keepmod~=0);

SRMt1keep=cell_registered_struct.cell_to_index_map(:,2);
SRMt1keepmod=SRMt1keep(SRMt2keep==0);
SRMt1keepmodSRMt1=SRMt1keepmod(SRMt1keepmod~=0);

t3keep=cell_registered_struct.cell_to_index_map(:,3);
t3keepmod=t3keep(SRMt2keep==0&SRMt1keep==0);
t3keepmodSRMt1=t3keepmod(t3keepmod~=0);

keeps={SRMt2keepmodSRMt1 SRMt1keepmodSRMt1 t3keepmodSRMt1};

neuronall.A=[];
neuronall.A_k=[];
neuronall.coor=[];
neuronall.coor_k=[];

neuronall.reg=cell_registered_struct;

for i=1:length(cell_registered_struct.spatial_footprints_corrected)
    curr=neuronall.reg.spatial_footprints_corrected{i};
    centroids=neuronall.reg.centroid_locations_corrected{i};
    currkeep=keeps{i};
    currsz=size(curr);
    temp=reshape(curr,currsz(1),currsz(2)*currsz(3))';
    temp_keep=temp(:,currkeep);
    neuronall.A=[neuronall.A temp];
    neuronall.A_k=[neuronall.A_k temp_keep];
    temp_coor=centroids(currkeep,:);
    neuronall.coor=[neuronall.coor; centroids];
    neuronall.coor_k=[neuronall.coor_k; temp_coor];
end



%% Get Vidframes
SRMt1=load('SRMt1');
SRMt2=load('SRMt2');
SRMt3=load('SRMt3');

SRMt1vidframes=SRMt1.vidframes;
SRMt2vidframes=SRMt2.vidframes;
t3vidframes=SRMt3.vidframes;

t1vfend=SRMt1vidframes{length(SRMt1vidframes)};
t1vfend=t1vfend(2);

t2vfend=SRMt2vidframes{length(SRMt2vidframes)};
t2vfend=t2vfend(2);

t3vfend=t3vidframes{length(t3vidframes)};
t3vfend=t3vfend(2);

vidframesall_bytrial={[1 t1vfend] [t1vfend+1 t2vfend+t1vfend] [t1vfend+t2vfend+1 t2vfend+t1vfend+t3vfend]};
vidframesall={SRMt1vidframes SRMt2vidframes t3vidframes};


%% Add in Dropped Frames

nname=neuronall;
numtests=3;

t1dropped={[]
};

t2dropped={
[]
};

t3dropped={
[]};


for i=1:length(t1dropped)
    curr=t1dropped{i};
    if ~isnan(curr)
        curr=curr+1;
    end
    t1dropped{i}=curr;
end

for i=1:length(t2dropped)
    curr=t2dropped{i};
    if ~isnan(curr)
        curr=curr+1;
    end
    t2dropped{i}=curr;
end

if numtests==3
for i=1:length(t3dropped)
    curr=t3dropped{i};
    if ~isnan(curr)
        curr=curr+1;
    end
    t3dropped{i}=curr;
end
end

n=0;
vidframesall_bytrial={};
for i=1:numtests
    currtri=vidframesall{i};
    currfirst=currtri{1};
    currlast=currtri{length(currtri)};
    start=currfirst(1);
    last=currlast(2);
    vidframesall_bytrial{i}=[n+start n+last];
    n=n+last;
end


dropped={t1dropped t2dropped t3dropped};


[neuronall.c_dropped,vfbt_dropped,vfall_dropped]=CorrectDroppedFrames(neuronall.C,dropped,numtests,vidframesall,vidframesall_bytrial);
[neuronall.yra_dropped,vfbt_dropped,vfall_dropped]=CorrectDroppedFrames(neuronall.YrA,dropped,numtests,vidframesall,vidframesall_bytrial);

%% Reorder Cells

t1nums=SRMt2keep>0;
t2nums=SRMt1keep>0&SRMt2keep==0;
t3nums=t3keep>0&SRMt1keep==0&SRMt2keep==0;

ids=1:length(t1nums);
t1ids=ids(t1nums==1);
t2ids=ids(t2nums==1);
t3ids=ids(t3nums==1);

neuronall_mod=neuronall;

neuronall_mod.c_dropped=[neuronall.c_dropped(t1ids,:); neuronall_mod.c_dropped(t2ids,:); neuronall_mod.c_dropped(t3ids,:)];
neuronall_mod.yra_dropped=[neuronall.yra_dropped(t1ids,:); neuronall_mod.yra_dropped(t2ids,:); neuronall_mod.yra_dropped(t3ids,:)];

%% Generate GUI objects

vidnum1=length(vidframesall{1});
vidnum2=length(vidframesall{2});
vidnum3=length(vidframesall{3});

neuron_SRMt1.all_A=sum(neuron_SRMt1.nw_A');
neuron_SRMt2.all_A=sum(neuron_SRMt2.nw_A');
neuron_SRMt3.all_A=sum(neuron_SRMt3.nw_A');

[nt1,nt2,nt3]=Format4GUI(neuron_SRMt1,neuron_SRMt2,neuron_SRMt3,cell_registered_struct,keeps,vidframesall,3,3,3);

%% Get Kept Neurons

nt1_mod=load('neuron_L1N1_out');
nt1_mod=nt1_mod.modneuron;
nt2_mod=load('neuron_LL_out');
nt2_mod=nt2_mod.modneuron;
nt3_mod=load('neuron_NN_out');
nt3_mod=nt3_mod.modneuron;

neunumt1=length(nt1.c(:,1));
neunumt2=length(nt2.c(:,1));
neunumt3=length(nt3.c(:,1));

nt1_mod.rejectnums=nt1_mod.rejectnums;
nt2_mod.rejectnums=nt2_mod.rejectnums+neunumt1;
nt3_mod.rejectnums=nt3_mod.rejectnums+neunumt1+neunumt2;

neuronall_out=neuronall_mod;
neuronall_out.rejnums=sort([nt1_mod.rejectnums nt2_mod.rejectnums nt3_mod.rejectnums]);

%% 

ids=1:length(neuronall_out.c_dropped(:,1));
keepnums=ids(~ismember(ids,neuronall_out.rejnums));
neuronall_adj.keepnums=keepnums;

neuronall_out.c_dropped=neuronall_mod.c_dropped(keepnums,:);
neuronall_out.yra_dropped=neuronall_mod.yra_dropped(keepnums,:);

%% Normalize Noise between sessions

a=vfall_dropped{1};
b=vfall_dropped{2};

neuronall_adj=neuronall_out;

for i=1:length(neuronall_out.c_dropped(:,1))
      if max(neuronall_out.c_dropped(i,vid1(1):vid1(2)))~=0
        t1yra=neuronall_out.yra_dropped(i,vid1(1):vid1(2));
        t1c=neuronall_out.c_dropped(i,vid1(1):vid1(2));
        t1c_noisy=t1c+t1yra;
        pophist=reshape(t1c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt1=std(pophist(silent==1));
    
        t2yra=neuronall_out.yra_dropped(i,vid2(1):vid2(2));
        t2c=neuronall_out.c_dropped(i,vid2(1):vid2(2));
        t2c_noisy=t2c+t2yra;
        pophist=reshape(t2c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt2=std(pophist(silent==1));
        
        adjt1t2=sigt1/sigt2;

        if isnan(adjt1t2)
            adjt1t2=1;
        end

        nwt2c_noisy=adjt1t2*t2c_noisy;
   
        t3yra=neuronall_out.yra_dropped(i,vid3(1):vid3(2));
        t3c=neuronall_out.c_dropped(i,vid3(1):vid3(2));
        t3c_noisy=t3c+t3yra;
        pophist=reshape(t3c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt3=std(pophist(silent==1));
        adjt1t3=sigt1/sigt3;
         
        if isnan(adjt1t3)
          adjt1t3=1;
        end
        
        nwt3c_noisy=adjt1t3*t3c_noisy;

        neuronall_adj.noisyC(i,:)=[t1c_noisy nwt2c_noisy nwt3c_noisy];
      else
        t1yra=neuronall_out.yra_dropped(i,vid1(1):vid1(2));
        t1c=neuronall_out.c_dropped(i,vid1(1):vid1(2));
        t1c_noisy=t1c+t1yra;
        t2yra=neuronall_out.yra_dropped(i,vid2(1):vid2(2));
        t2c=neuronall_out.c_dropped(i,vid2(1):vid2(2));
        t2c_noisy=t2c+t2yra;
        
        pophist=reshape(t2c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt2=std(pophist(silent==1));
        
        t3yra=neuronall_out.yra_dropped(i,vid3(1):vid3(2));
        t3c=neuronall_out.c_dropped(i,vid3(1):vid3(2));
        t3c_noisy=t3c+t3yra;
        pophist=reshape(t3c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt3=std(pophist(silent==1));
        
        adjt2t3=sigt2/sigt3;
        
        if isnan(adjt2t3)
            adjt2t3=1;
        end
    
        nwt3c_noisy=t3c_noisy*adjt2t3;
        neuronall_adj.noisyC(i,:)=[t1c_noisy t2c_noisy nwt3c_noisy];
   
      end
end

%% Set constant baseline 
neuronall_adj.adjnoisyC=neuronall_adj.noisyC;

fr={vid1 vid2 vid3};
timebin=1000;
   
for i=1:length(neuronall_adj.adjnoisyC(:,1))
    baseline=0;
    for sess=1:length(fr)
    currfr=fr{sess};
    currtri=neuronall_adj.adjnoisyC(i,currfr(1):currfr(2));
    numbins=floor(length(neuronall_adj.adjnoisyC(i,currfr(1):currfr(2)))/timebin);
    for bin=1:numbins
        if bin<numbins
            currbin=neuronall_adj.adjnoisyC(i,(currfr(1)+timebin*(bin-1)):(currfr(1)-1+timebin*bin));
            currmed=median(currbin(~isnan(currbin)));
            nwbin=currbin-(currmed-baseline);
            neuronall_adj.adjnoisyC_v2(i,(currfr(1)+timebin*(bin-1)):(currfr(1)-1+timebin*bin))=nwbin;
        else
            currbin=neuronall_adj.adjnoisyC(i,(currfr(1)+timebin*(bin-1)):currfr(2));
            currmed=median(currbin(~isnan(currbin)));
            nwbin=currbin-(currmed-baseline);
            neuronall_adj.adjnoisyC_v2(i,(currfr(1)+timebin*(bin-1)):currfr(2))=nwbin;
        end
    end
        
    end
end
%% adjust noise again

a=vfall_dropped{1};
b=vfall_dropped{2};

for i=1:length(neuronall_adj.adjnoisyC.c_dropped(:,1))
      if max(neuronall_adj.adjnoisyC.c_dropped(i,vid1(1):vid1(2)))~=0
        t1yra=neuronall_adj.adjnoisyC.yra_dropped(i,vid1(1):vid1(2));
        t1c=neuronall_adj.adjnoisyC.c_dropped(i,vid1(1):vid1(2));
        t1c_noisy=t1c+t1yra;
        pophist=reshape(t1c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt1=std(pophist(silent==1));
    
        t2yra=neuronall_adj.adjnoisyC.yra_dropped(i,vid2(1):vid2(2));
        t2c=neuronall_adj.adjnoisyC.c_dropped(i,vid2(1):vid2(2));
        t2c_noisy=t2c+t2yra;
        pophist=reshape(t2c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt2=std(pophist(silent==1));
        
        adjt1t2=sigt1/sigt2;

        if isnan(adjt1t2)
            adjt1t2=1;
        end

        nwt2c_noisy=adjt1t2*t2c_noisy;
   
        t3yra=neuronall_adj.adjnoisyC.yra_dropped(i,vid3(1):vid3(2));
        t3c=neuronall_adj.adjnoisyC.c_dropped(i,vid3(1):vid3(2));
        t3c_noisy=t3c+t3yra;
        pophist=reshape(t3c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt3=std(pophist(silent==1));
        adjt1t3=sigt1/sigt3;
         
        if isnan(adjt1t3)
          adjt1t3=1;
        end
        
        nwt3c_noisy=adjt1t3*t3c_noisy;

        neuronall_adj.adjnoisyC_v2(i,:)=[t1c_noisy nwt2c_noisy nwt3c_noisy];
      else
        t1yra=neuronall_adj.adjnoisyC.yra_dropped(i,vid1(1):vid1(2));
        t1c=neuronall_adj.adjnoisyC.c_dropped(i,vid1(1):vid1(2));
        t1c_noisy=t1c+t1yra;
        t2yra=neuronall_adj.adjnoisyC.yra_dropped(i,vid2(1):vid2(2));
        t2c=neuronall_adj.adjnoisyC.c_dropped(i,vid2(1):vid2(2));
        t2c_noisy=t2c+t2yra;
        
        pophist=reshape(t2c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt2=std(pophist(silent==1));
        
        t3yra=neuronall_adj.adjnoisyC.yra_dropped(i,vid3(1):vid3(2));
        t3c=neuronall_adj.adjnoisyC.c_dropped(i,vid3(1):vid3(2));
        t3c_noisy=t3c+t3yra;
        pophist=reshape(t3c_noisy,[],1);
        pophist=pophist(~isnan(pophist));
        pop_offset=quantile(pophist,[0.05 0.95]);
        silent=pophist<pop_offset(2)&pophist>pop_offset(1);
        sigt3=std(pophist(silent==1));
        
        adjt2t3=sigt2/sigt3;
        
        if isnan(adjt2t3)
            adjt2t3=1;
        end
    
        nwt3c_noisy=t3c_noisy*adjt2t3;
        neuronall_adj.adjnoisyC_v2.noisyC(i,:)=[t1c_noisy t2c_noisy nwt3c_noisy];
   
      end
end

%% Set constant baseline 2
neuronall_adj.adjnoisyC_v3=neuronall_adj.adjnoisyC_v2;
fr={vid1 vid2 vid3};

for i=1:length(neuronall_adj.adjnoisyC_v3(:,1))
    baseline=0;
    for sess=1:length(fr)
    currfr=fr{sess};
    currtri=neuronall_adj.adjnoisyC_v3(i,currfr(1):currfr(2));
    if max(currtri)~=0
    numbins=floor(length(neuronall_adj.adjnoisyC_v3(i,currfr(1):currfr(2)))/1000);
    for bin=1:numbins
        if bin<numbins
            currbin=neuronall_adj.adjnoisyC_v3(i,(currfr(1)+1000*(bin-1)):(currfr(1)-1+1000*bin));
            currmed=median(currbin(~isnan(currbin)));
            nwbin=currbin-(currmed-baseline);
            neuronall_adj.adjnoisyC_v3(i,(currfr(1)+1000*(bin-1)):(currfr(1)-1+1000*bin))=nwbin;
        else
            currbin=neuronall_adj.adjnoisyC_v3(i,(currfr(1)+1000*(bin-1)):currfr(2));
            currmed=median(currbin(~isnan(currbin)));
            nwbin=currbin-(currmed-baseline);
            neuronall_adj.adjnoisyC_v3(i,(currfr(1)+1000*(bin-1)):currfr(2))=nwbin;
        end
    end       
    end
end
        

%% Get Deconvolved signal
neuronall_decon=neuronall_adj;
neuronall_decon.noisyC=neuronall_adj.adjnoisyC;
neuronall_decon.noisyC(isnan(neuronall_adj.noisyC))=0;

vid1=vfbt_dropped{1};
vid2=vfbt_dropped{2};
vid3=vfbt_dropped{3};


for i=1:length(neuronall_decon.noisyC(:,1))
    t1c=neuronall_decon.noisyC(i,vid1(1):vid1(2));
    t2c=neuronall_decon.noisyC(i,vid2(1):vid2(2));
    t3c=neuronall_decon.noisyC(i,vid3(1):vid3(2));
    
    currctrace=[];
    currstrace=[];
    
    if max(t1c)~=0&&max(t2c)~=0&&max(t3c)~=0
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        [ccurr,scurr,options]=deconvolveCa(t1c,'thresholded','ar2','b',median(t1c),'optimize_sn','TRUE','smin',mad(t1c,1),'maxIter',100,'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        [ccurr,scurr]=deconvolveCa(t2c,'thresholded','ar2','b',median(t2c),'optimize_sn','TRUE','smin',mad(t2c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        [ccurr,scurr]=deconvolveCa(t3c,'thresholded','ar2','b',median(t3c),'optimize_sn','TRUE','smin',mad(t3c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
    elseif max(t1c)==0&&max(t2c)~=0&&max(t3c)~=0
        currctrace=t1c;
        currstrace=t1c;
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        [ccurr,scurr]=deconvolveCa(t2c,'thresholded','ar2','b',median(t2c),'optimize_sn','TRUE','smin',mad(t2c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        [ccurr,scurr]=deconvolveCa(t3c,'thresholded','ar2','b',median(t3c),'optimize_sn','TRUE','smin',mad(t3c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
    elseif max(t1c)==0&&max(t2c)==0&&max(t3c)~=0
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        currctrace=[t1c t2c];
        currstrace=[t1c t2c];
        [ccurr,scurr]=deconvolveCa(t3c,'thresholded','ar2','b',median(t3c),'optimize_sn','TRUE','smin',mad(t3c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
    elseif max(t1c)~=0&&max(t2c)==0&&max(t3c)==0
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        [ccurr,scurr,options]=deconvolveCa(t1c,'thresholded','ar2','b',median(t1c),'optimize_sn','TRUE','smin',mad(t1c,1),'maxIter',100,'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        currctrace=[currctrace t2c];
        currstrace=[currstrace t2c];
        currctrace=[currctrace t3c];
        currstrace=[currstrace t3c];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
    elseif max(t1c)~=0&&max(t2c)~=0&&max(t3c)==0
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        [ccurr,scurr,options]=deconvolveCa(t1c,'thresholded','ar2','b',median(t1c),'optimize_sn','TRUE','smin',mad(t1c,1),'maxIter',100,'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        [ccurr,scurr]=deconvolveCa(t2c,'thresholded','ar2','b',median(t2c),'optimize_sn','TRUE','smin',mad(t2c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        currctrace=[currctrace t3c];
        currstrace=[currstrace t3c];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
     elseif max(t1c)~=0&&max(t2c)==0&&max(t3c)~=0
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        [ccurr,scurr,options]=deconvolveCa(t1c,'thresholded','ar2','b',median(t1c),'optimize_sn','TRUE','smin',mad(t1c,1),'maxIter',100,'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        currctrace=[currctrace t2c];
        currstrace=[currstrace t2c];
        [ccurr,scurr]=deconvolveCa(t3c,'thresholded','ar2','b',median(t3c),'optimize_sn','TRUE','smin',mad(t3c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
     elseif max(t1c)==0&&max(t2c)~=0&&max(t3c)==0
        trivf1=vfall_dropped{1};
        trivf2=vfall_dropped{2};
        trivf3=vfall_dropped{3};
        currctrace=t1c;
        currstrace=t1c;
        [ccurr,scurr]=deconvolveCa(t2c,'thresholded','ar2','b',median(t2c),'optimize_sn','TRUE','smin',mad(t2c,1),'maxIter',100, 'thresh_factor',2);
        currctrace=[currctrace ccurr'];
        currstrace=[currstrace scurr'];
        currctrace=[currctrace t3c];
        currstrace=[currstrace t3c];
        neuronall_decon.c(i,:)=currctrace;
        neuronall_decon.s(i,:)=currstrace;
    end
end


%%
neuronall_decon.coor_k=neuronall_out.coor_k(keepnums,:);
neuronall_decon.A_k=neuronall_out.A_k(:,keepnums);
neuronall_decon.all_A=sum(neuronall_decon.A_k');
neuronall_decon.dims=[adjusted_y_size adjusted_x_size];
%% Clear previous

%You will want to export your neuron_decon object, as well as:
%vfall_dropped, vfbt_dropped, dims, and droppedframes

clear all;

%% Align with behavior

filename='neuron_t1t2t3';
data=load(filename);

neuron=data.neuronall_decon;
vidframesall=data.vfall_dropped;
vidframes_bytrial=data.vfbt_dropped;
dropped=data.dropped;

SRMt1=data.SRMt1;
SRMt2=data.SRMt2;
SRMt3=data.SRMt3;

%% Align neural data and behavioral data

numvidt1=3;
numvidt2=3;
numvidt3=3;

trilng=[];
prompt='Input original trial lengths';
for i=1:(numvidt1+numvidt2+numvidt3)
    trilngtemp=input(prompt,'s');
    trilngtemp=str2double(trilngtemp);
    trilng(i)=trilngtemp;
end

modnums=[];
prompt='Input modnums if any';
for i=1:(numvidt1+numvidt2+numvidt3)
    modtemp=input(prompt,'s');
    modtemp=str2double(modtemp);
    modnums(i)=modtemp;
end

lims=[];
prompt='Input lims, if any';
for i=1:(numvidt1+numvidt2+numvidt3)
    limstemp=input(prompt,'s');
    limstemp=str2double(limstemp);
    lims(i)=limstemp;
end

skipt1=[];
skipt2=[];
skipt3=[];
prompt='Input which test 1 trials, if any, to skip';
skipt1=input(prompt,'s');
skipt1=str2num(skipt1);
prompt='Input which test 2 trials, if any, to skip';
skipt2=input(prompt,'s');
skipt2=str2num(skipt2);
if t3==1
prompt='Input which test 3 trials, if any, to skip';
skipt3=input(prompt,'s');
skipt3=str2num(skipt3);
end


%% Put in neut1

triframesall={SRMt1.triframes SRMt2.triframes SRMt3.triframes};

t1vfend=SRMt1.triframes{length(SRMt1.triframes)};
t1vfend=t1vfend(2);

t2vfend=SRMt2.triframes{length(SRMt2.triframes)};
t2vfend=t2vfend(2);

t3vfend=SRMt3.triframes{length(SRMt3.triframes)};
t3vfend=t3vfend(2);
tf3=[1 t3vfend];
triframes_bytrial={[1 t1vfend] [t1vfend+1 t2vfend+t1vfend] [t1vfend+t2vfend+1 t2vfend+t1vfend+t3vfend]};
tfall3=triframesall{3};

currtf=[];

tf1=[1 t1vfend];
tf2=[1 t2vfend];
tf3=[1 t3vfend];

tfall1=triframesall{1};
tfall2=triframesall{2};
tfall3=triframesall{3};

vf1all=vidframesall{1};
vf2all=vidframesall{2};
vf3all=vidframesall{3};

vfbt1=vidframes_bytrial{1};
neut1.c=neuron.c(:,vfbt1(1):vfbt1(2));
neut1.c_beh=[];
neut1.c_raw=neuron.adjnoisyC_v4(:,vfbt1(1):vfbt1(2));
neut1.craw_beh=[];
neut1.s=neuron.s(:,vfbt1(1):vfbt1(2));
neut1.s_beh=[];

neut1.adjDLxytrack=[];
neut1.adjtriframes=SRMt1.SRMt1.triframes;
neut1.adjDLxytrackb=[];
neut1.adjDLcircbeh=[];
neut1.adjbeh=[];
neut1.adjanybeh=[];
neut1.adjxytrack=[];

vfbt2=vidframes_bytrial{2};
neut2.c=neuron.c(:,vfbt2(1):vfbt2(2));
neut2.c_beh=[];
neut2.c_raw=neuron.adjnoisyC_v4(:,vfbt2(1):vfbt2(2));
neut2.craw_beh=[];
neut2.s=neuron.s(:,vfbt2(1):vfbt2(2));
neut2.s_beh=[];

if t3==1
    vfbt3=vidframes_bytrial{3};
    neut3.c=neuron.c(:,vfbt3(1):vfbt3(2));
    neut3.c_beh=[];
    neut3.c_raw=neuron.adjnoisyC_v4(:,vfbt3(1):vfbt3(2));
    neut3.craw_beh=[];
    neut3.s=neuron.s(:,vfbt3(1):vfbt3(2));
    neut3.s_beh=[];
end


%% Replace Dropped with NaN;

t1dropped=dropped{1};
t2dropped=dropped{2};
t3dropped=dropped{3};
vfbt1=vidframes_bytrial{1};
vfbt2=vidframes_bytrial{2};
vfbt3=vidframes_bytrial{3};


for i=1:length(t1dropped)
    if ~isnan(t1dropped{i})
        cf=vf1all{i};
        currneu.c=neut1.c(:,cf(1):cf(2));
        currneu.s=neut1.s(:,cf(1):cf(2));
        currneu.c_raw=neut1.c_raw(:,cf(1):cf(2));
        currdropped=t1dropped{i};
        for d=1:length(currdropped)
            currneu.c(:,currdropped(d))=NaN(length(currneu.c(:,1)),1);
            currneu.s(:,currdropped(d))=NaN(length(currneu.s(:,1)),1);
            currneu.c_raw(:,currdropped(d))=NaN(length(currneu.c_raw(:,1)),1);
        end
        neut1.c(:,cf(1):cf(2))=currneu.c;
        neut1.c_raw(:,cf(1):cf(2))=currneu.c_raw;
        neut1.s(:,cf(1):cf(2))=currneu.s;
    end
end

for i=1:length(t2dropped)
    if ~isnan(t2dropped{i})
        cf=vf2all{i};
        currneu.c=neut2.c(:,cf(1):cf(2));
        currneu.s=neut2.s(:,cf(1):cf(2));
        currneu.c_raw=neut2.c_raw(:,cf(1):cf(2));
        currdropped=t2dropped{i};
        for d=1:length(currdropped)
            currneu.c(:,currdropped(d))=NaN(length(currneu.c(:,1)),1);
            currneu.s(:,currdropped(d))=NaN(length(currneu.s(:,1)),1);
            currneu.c_raw(:,currdropped(d))=NaN(length(currneu.c_raw(:,1)),1);
        end
        neut2.c(:,cf(1):cf(2))=currneu.c;
        neut2.c_raw(:,cf(1):cf(2))=currneu.c_raw;
        neut2.s(:,cf(1):cf(2))=currneu.s;
    end
end

for i=1:length(t3dropped)
    if ~isnan(t3dropped{i})
        cf=vf3all{i};
        currneu.c=neut3.c(:,cf(1):cf(2));
        currneu.s=neut3.s(:,cf(1):cf(2));
        currneu.c_raw=neut3.c_raw(:,cf(1):cf(2));
        currdropped=t3dropped{i};
        for d=1:length(currdropped)
            currneu.c(:,currdropped(d))=NaN(length(currneu.c(:,1)),1);
            currneu.s(:,currdropped(d))=NaN(length(currneu.s(:,1)),1);
            currneu.c_raw(:,currdropped(d))=NaN(length(currneu.c_raw(:,1)),1);
        end
        neut3.c(:,cf(1):cf(2))=currneu.c;
        neut3.c_raw(:,cf(1):cf(2))=currneu.c_raw;
        neut3.s(:,cf(1):cf(2))=currneu.s;
    end
end

neuron.s(:,vfbt1(1):vfbt1(2))=neut1.s;
neuron.c(:,vfbt1(1):vfbt1(2))=neut1.c;
neuron.c_raw(:,vfbt1(1):vfbt1(2))=neut1.c_raw;

neuron.s(:,vfbt2(1):vfbt2(2))=neut2.s;
neuron.c(:,vfbt2(1):vfbt2(2))=neut2.c;
neuron.c_raw(:,vfbt2(1):vfbt2(2))=neut2.c_raw;

neuron.s(:,vfbt3(1):vfbt3(2))=neut3.s;
neuron.c(:,vfbt3(1):vfbt3(2))=neut3.c;
neuron.c_raw(:,vfbt3(1):vfbt3(2))=neut3.c_raw;


%% 
numvidt1=length(vf1all);
numvidt2=length(vf2all);

neuron.craw_beh=[];
neuron.c_beh=[];
neuron.s_beh=[];

numvidt3=length(vf3all);
nwtf3=tfall3;

nwtf1=tfall1;
nwtf2=tfall2;

tri=1;
currnum=1;
start=1;

for firstt=1:numvidt1   
    if ~ismember(firstt,skipt1)
        currtf=tfall1{tri};
        tri=tri+1;
        currvf=vf1all{firstt};
        
        if isnan(mean(lims{currnum}))
            if ~isnan(modnums(currnum))
                modinit=modnums(currnum);
           numvidt1=length(vf1all);
            elseif (currvf(2)-currvf(1))>trilng(currnum)   %If the calcium imaging video extended beyond the behavior video, find end of trial
                modinit=currvf(2)-currvf(1)+1-trilng(currnum);
            else
                modinit=0;
            end
            if (currvf(2)-currvf(1))>(currtf(2)-currtf(1))
                modvid=[(currvf(2)-(currtf(2)-currtf(1))-modinit) (currvf(2)-modinit)]; %behcups1(2)-behcups1(1) = Length of behavior trial; From end of cups calcium imaging video, subtract the length of the behavior trial and any extra frames
                modbvid=currtf;
            else
                modvid=currvf;
                modbvid=[(currtf(2)-(currvf(2)-currvf(1))-modinit) (currtf(2)-modinit)]; %behcups1(2)-behcups1(1) = Length of behavior trial; From end of cups calcium imaging video, subtract the length of the behavior trial and any extra frames
            end
        else
            prevvf=vf1all{tri-1};
            currlims=lims{currnum};
            modvid(1)=prevvf(2)+1+currlims(1);
            modvid(2)=modvid(1)+currlims(2);
            nwtf1{firstt}=[currtf(1) currtf(1)+currlims(2)];
            for rem=(firstt+1):numvidt1
                frames2mod=nwtf1{rem};
                diff=currtf(2)-currtf(1)-currlims(2);
                frames2mod=frames2mod-diff;
                nwtf1{rem}=frames2mod;
            end 
        end
        
            neut1.adjbeh=[neut1.adjbeh SRMt1.SRMt1.beh(:,modbvid(1):modbvid(2))];
            neut1.adjanybeh=[neut1.adjanybeh SRMt1.SRMt1.anybeh(:,modbvid(1):modbvid(2))];
            neut1.adjxytrack=[neut1.adjxytrack SRMt1.SRMt1.xytrack(:,modbvid(1):modbvid(2))];;
            neut1.adjDLxytrack=[neut1.adjDLxytrack SRMt1.SRMt1.DLxytrack(:,modbvid(1):modbvid(2))];;
            neut1.adjDLxytrackb=[neut1.adjDLxytrackb SRMt1.SRMt1.DLxytrackb(:,modbvid(1):modbvid(2))];;
            neut1.adjDLcircbeh=[neut1.adjDLcircbeh SRMt1.SRMt1.DLcircbeh(:,modbvid(1):modbvid(2))];
            currfr=SRMt1.SRMt1.triframes{firstt};
            currlng=modbvid(2)-modbvid(1);
            last=start+currlng;
            nwfr=[start last];
            start=last+1;
            neut1.adjtriframes{firstt}=nwfr;

        neut1.c_beh=[neut1.c_beh neut1.c(:,modvid(1):modvid(2))];
        neut1.s_beh=[neut1.s_beh neut1.s(:,modvid(1):modvid(2))];
        neut1.craw_beh=[neut1.craw_beh neut1.c_raw(:,modvid(1):modvid(2))];
        currnum=currnum+1;
    end
end


    
        
tri=1;
for secondt=1:numvidt2
    if ~ismember(secondt,skipt2)
        currtf=tfall2{tri};
        tri=tri+1;
        currvf=vf2all{secondt};
        
        if isnan(mean(lims{currnum}))
            if ~isnan(modnums(currnum))
                modinit=modnums(currnum);
            elseif (currvf(2)-currvf(1))>trilng(currnum)   %If the calcium imaging video extended beyond the behavior video, find end of trial
                modinit=currvf(2)-currvf(1)-trilng(currnum);
            else
                modinit=0;
            end
        
            modvid=[(currvf(2)-(currtf(2)-currtf(1))-modinit) (currvf(2)-modinit)]; %behcups1(2)-behcups1(1) = Length of behavior trial; From end of cups calcium imaging video, subtract the length of the behavior trial and any extra frames
        else
            prevvf=vf2all{tri-1};
            currlims=lims{currnum};
            modvid(1)=prevvf(2)+1+currlims(1);
            modvid(2)=modvid(1)+currlims(2);
            nwtf2{secondt}=[currtf(1) currtf(1)+currlims(2)];
            for rem=secondt:numvidt2
                frames2mod=nwtf2{rem};
                diff=currtf(2)-currtf(1)-currlims(2);
                frames2mod=frames2mod-diff;
                nwtf2{rem}=frames2mod;
            end 
        end
        
        neut2.c_beh=[neut2.c_beh neut2.c(:,modvid(1):modvid(2))];
        neut2.s_beh=[neut2.s_beh neut2.s(:,modvid(1):modvid(2))];
        neut2.craw_beh=[neut2.craw_beh neut2.c_raw(:,modvid(1):modvid(2))];
        currnum=currnum+1;
    end
        
end

tri=1;
for thirdt=1:numvidt3
    if ~ismember(thirdt,skipt3)
        currtf=tfall3{tri};
        tri=tri+1;
        currvf=vf3all{thirdt};
        
        if isnan(mean(lims{currnum}))
            if ~isnan(modnums(currnum))
                modinit=modnums(currnum);
            elseif (currvf(2)-currvf(1))>trilng(currnum)   %If the calcium imaging video extended beyond the behavior video, find end of trial
                modinit=currvf(2)-currvf(1)-trilng(currnum);
            else
                modinit=0;
            end
            modvid=[(currvf(2)-(currtf(2)-currtf(1))-modinit) (currvf(2)-modinit)]; %behcups1(2)-behcups1(1) = Length of behavior trial; From end of cups calcium imaging video, subtract the length of the behavior trial and any extra frames
        else
            prevvf=vf3all{tri-1};
            currlims=lims{currnum};
            modvid(1)=prevvf(2)+1+currlims(1);
            modvid(2)=modvid(1)+currlims(2);
            nwtf3{thirdt}=[currtf(1) currtf(1)+currlims(2)];
            for rem=thirdt:numvidt3
                frames2mod=nwtf3{rem};
                diff=currtf(2)-currtf(1)-currlims(2);
                frames2mod=frames2mod-diff;
                nwtf3{rem}=frames2mod;
            end 
        end
        
        neut3.c_beh=[neut3.c_beh neut3.c(:,modvid(1):modvid(2))];
        neut3.s_beh=[neut3.s_beh neut3.s(:,modvid(1):modvid(2))];
        neut3.craw_beh=[neut3.craw_beh neut3.c_raw(:,modvid(1):modvid(2))];
       currnum=currnum+1;
    end
end


neuron.c_beh=[neut1.c_beh neut2.c_beh neut3.c_beh];
neuron.s_beh=[neut1.s_beh neut2.s_beh neut3.s_beh];
neuron.craw_beh=[neut1.craw_beh neut2.craw_beh neut3.craw_beh];



%% Combinetriframes

neuron.triframes=nwtf1;

numtrit1=3;
numtrit2=3;
if t3==1
    numtrit3=3;
end
if t4==1
    numtrit4=3;
    numtrit5=3;
end

lasttrit1=nwtf1{numtrit1};
lastt1=lasttrit1(2);

lasttrit2=nwtf2{numtrit2};
lastt2=lasttrit2(2);

lasttrit3=nwtf3{numtrit3};
lastt3=lasttrit3(2);

DLxytrackall=[];
DLxytrackallb=[];
DLcircbehall=[];
behall=[];

for i=1:numtrit1
    currfr=neuron.triframes{i};
    lng=currfr(2)-currfr(1);
    orgfr=SRMt1.SRMt1.triframes{i};
    DLxytrackall=[DLxytrackall SRMt1.SRMt1.DLxytrack(:,orgfr(1):orgfr(1)+lng)];
    DLxytrackallb=[DLxytrackallb SRMt1.SRMt1.DLxytrackb(:,orgfr(1):orgfr(1)+lng)];
    DLcircbehall=[DLcircbehall SRMt1.SRMt1.DLcircbeh(:,orgfr(1):orgfr(1)+lng)];
    behall=[behall SRMt1.SRMt1.beh(:,orgfr(1):orgfr(1)+lng)];
end

for i=1:numtrit2
    currfr=nwtf2{i};
    lng=currfr(2)-currfr(1);
    orgfr=SRMt2.SRMt2.triframes{i};
    DLxytrackall=[DLxytrackall SRMt2.SRMt2.DLxytrack(:,orgfr(1):orgfr(1)+lng)];
    DLxytrackallb=[DLxytrackallb SRMt2.SRMt2.DLxytrackb(:,orgfr(1):orgfr(1)+lng)];
    DLcircbehall=[DLcircbehall SRMt2.SRMt2.DLcircbeh(:,orgfr(1):orgfr(1)+lng)];
    behall=[behall SRMt2.SRMt2.beh(:,orgfr(1):orgfr(1)+lng)];
    currfr=currfr+lastt1;
    neuron.triframes{i+numtrit1}=currfr;
end

for i=1:length(nwtf3)
    currfr=nwtf3{i};
    lng=currfr(2)-currfr(1);
    orgfr=SRMt3.SRMt3.triframes{i};
    DLxytrackall=[DLxytrackall SRMt3.SRMt3.DLxytrack(:,orgfr(1):orgfr(1)+lng)];
    DLxytrackallb=[DLxytrackallb SRMt3.SRMt3.DLxytrackb(:,orgfr(1):orgfr(1)+lng)];
    DLcircbehall=[DLcircbehall SRMt3.SRMt3.DLcircbeh(:,orgfr(1):orgfr(1)+lng)];
    behall=[behall SRMt3.SRMt3.beh(:,orgfr(1):orgfr(1)+lng)];
    currfr=currfr+lastt1+lastt2;
    neuron.triframes{i+numtrit1+numtrit2}=currfr;
end


t=DLcircbehall(1,:);
DLcircbehall=DLcircbehall(2:4,:);
behall=behall(2:3,:);
%% Save Data to H5

h5create('msID_SRMt1.h5','/neuron/C',size(neuronSRMt1.c_beh'));
h5create('msID_SRMt1.h5','/neuron/S',size(neuronSRMt1.s_beh'));
h5create('msID_SRMt2.h5','/neuron/C',size(neuronSRMt2.c_beh'));
h5create('msID_SRMt2.h5','/neuron/S',size(neuronSRMt2.s_beh'));
h5create('msID_SRMt3.h5','/neuron/C',size(neuronSRMt3.c_beh'));
h5create('msID_SRMt3.h5','/neuron/S',size(neuronSRMt3.s_beh'));
h5create('msID_SRMt1.h5','/neuron/C_Raw',size(neuronSRMt1.craw_beh'));
h5create('msID_SRMt2.h5','/neuron/C_Raw',size(neuronSRMt2.craw_beh'));
h5create('msID_SRMt3.h5','/neuron/C_Raw',size(neuronSRMt3.craw_beh'));
h5create('msID_SRMt1.h5','/neuron/beh/DLcircbeh',size(neuronSRMt1.DLcircbeh'));
h5create('msID_SRMt1.h5','/neuron/beh/scoredbeh',size(neuronSRMt1.beh'));
h5create('msID_SRMt1.h5','/neuron/MotionTrack/head',size(neuronSRMt1.DLxytrack'));
h5create('msID_SRMt1.h5','/neuron/MotionTrack/body',size(neuronSRMt1.DLxytrackb'));
h5create('msID_SRMt2.h5','/neuron/beh/DLcircbeh',size(neuronSRMt2.DLcircbeh'));
h5create('msID_SRMt2.h5','/neuron/beh/scoredbeh',size(neuronSRMt2.beh'));
h5create('msID_SRMt2.h5','/neuron/MotionTrack/head',size(neuronSRMt2.DLxytrack'));
h5create('msID_SRMt2.h5','/neuron/MotionTrack/body',size(neuronSRMt2.DLxytrackb'));
h5create('msID_SRMt3.h5','/neuron/beh/DLcircbeh',size(neuronSRMt3.DLcircbeh'));
h5create('msID_SRMt3.h5','/neuron/beh/scoredbeh',size(neuronSRMt3.beh'));
h5create('msID_SRMt3.h5','/neuron/MotionTrack/head',size(neuronSRMt3.DLxytrack'));
h5create('msID_SRMt3.h5','/neuron/MotionTrack/body',size(neuronSRMt3.DLxytrackb'));
h5create('msID_SRMt1.h5','/numtrials',size(numtrit1));
h5create('msID_SRMt2.h5','/numtrials',size(numtrit2));
h5create('msID_SRMt3.h5','/numtrials',size(numtrit3));
h5create('msID_SRMt1.h5','/ms1posn',size(SRMt1.SRMt1.L1posns));
h5create('msID_SRMt2.h5','/ms1posn',size(SRMt2.SRMt2.L1posns));
h5create('msID_SRMt3.h5','/ms1posn',size(SRMt3.SRMt3.L1posns));
h5create('msID_SRMt1.h5','/frames/Test1cups',size(neuronSRMt1.triframes{1}));
h5create('msID_SRMt1.h5','/frames/Test1novelvsfamiliar',size(neuronSRMt1.triframes{2}));
h5create('msID_SRMt1.h5','/frames/Test1novelvsfamiliar_posnswap',size(neuronSRMt1.triframes{3}));
h5create('msID_SRMt2.h5','/frames/Test2cups',size(neuronSRMt2.triframes{1}));
h5create('msID_SRMt2.h5','/frames/Test2novels',size(neuronSRMt2.triframes{2}));
h5create('msID_SRMt2.h5','/frames/Test2novels_posnswap',size(neuronSRMt2.triframes{3}));
h5create('msID_SRMt3.h5','/frames/Test3cups',size(neuronSRMt3.triframes{1}));
h5create('msID_SRMt3.h5','/frames/Test3littermates',size(neuronSRMt3.triframes{2}));
h5create('msID_SRMt3.h5','/frames/Test3littermates_posnswap',size(neuronSRMt3.triframes{3}));
h5create('msID_SRMt1.h5','/neuron/t',size(neuronSRMt1.t));
h5create('msID_SRMt2.h5','/neuron/t',size(neuronSRMt2.t));
h5create('msID_SRMt3.h5','/neuron/t',size(neuronSRMt3.t));
h5create('msID_SRMt1.h5','/neuron/IDarray',size(neuronSRMt1.IDarray),'DataType','string');
h5create('msID_SRMt2.h5','/neuron/IDarray',size(neuronSRMt2.IDarray),'DataType','string');
h5create('msID_SRMt3.h5','/neuron/IDarray',size(neuronSRMt3.IDarray),'DataType','string');

h5write('msID_SRMt1.h5','/neuron/C',neuronSRMt1.c_beh');
h5write('msID_SRMt1.h5','/neuron/S',neuronSRMt1.s_beh');
h5write('msID_SRMt2.h5','/neuron/C',neuronSRMt2.c_beh');
h5write('msID_SRMt2.h5','/neuron/S',neuronSRMt2.s_beh');
h5write('msID_SRMt3.h5','/neuron/C',neuronSRMt3.c_beh');
h5write('msID_SRMt3.h5','/neuron/S',neuronSRMt3.s_beh');
h5write('msID_SRMt1.h5','/neuron/C_Raw',neuronSRMt1.craw_beh');
h5write('msID_SRMt2.h5','/neuron/C_Raw',neuronSRMt2.craw_beh');
h5write('msID_SRMt3.h5','/neuron/C_Raw',neuronSRMt3.craw_beh');
h5write('msID_SRMt1.h5','/neuron/beh/DLcircbeh',neuronSRMt1.DLcircbeh');
h5write('msID_SRMt1.h5','/neuron/beh/scoredbeh',neuronSRMt1.beh');
h5write('msID_SRMt1.h5','/neuron/MotionTrack/head',neuronSRMt1.DLxytrack');
h5write('msID_SRMt1.h5','/neuron/MotionTrack/body',neuronSRMt1.DLxytrackb');
h5write('msID_SRMt2.h5','/neuron/beh/DLcircbeh',neuronSRMt2.DLcircbeh');
h5write('msID_SRMt2.h5','/neuron/beh/scoredbeh',neuronSRMt2.beh');
h5write('msID_SRMt2.h5','/neuron/MotionTrack/head',neuronSRMt2.DLxytrack');
h5write('msID_SRMt2.h5','/neuron/MotionTrack/body',neuronSRMt2.DLxytrackb');
h5write('msID_SRMt3.h5','/neuron/beh/DLcircbeh',neuronSRMt3.DLcircbeh');
h5write('msID_SRMt3.h5','/neuron/beh/scoredbeh',neuronSRMt3.beh');
h5write('msID_SRMt3.h5','/neuron/MotionTrack/head',neuronSRMt3.DLxytrack');
h5write('msID_SRMt3.h5','/neuron/MotionTrack/body',neuronSRMt3.DLxytrackb');
h5write('msID_SRMt1.h5','/numtrials',numtrit1);
h5write('msID_SRMt2.h5','/numtrials',numtrit2);
h5write('msID_SRMt3.h5','/numtrials',numtrit3);
h5write('msID_SRMt1.h5','/ms1posn',SRMt1.SRMt1.L1posns);
h5write('msID_SRMt2.h5','/ms1posn',SRMt2.SRMt2.L1posns);
h5write('msID_SRMt3.h5','/ms1posn',SRMt3.SRMt3.L1posns);
h5write('msID_SRMt1.h5','/frames/Test1cups',neuronSRMt1.triframes{1});
h5write('msID_SRMt1.h5','/frames/Test1novelvsfamiliar',neuronSRMt1.triframes{2});
h5write('msID_SRMt1.h5','/frames/Test1novelvsfamiliar_posnswap',neuronSRMt1.triframes{3});
h5write('msID_SRMt2.h5','/frames/Test2cups',neuronSRMt2.triframes{1});
h5write('msID_SRMt2.h5','/frames/Test2novels',neuronSRMt2.triframes{2});
h5write('msID_SRMt2.h5','/frames/Test2novels_posnswap',neuronSRMt2.triframes{3});
h5write('msID_SRMt3.h5','/frames/Test3cups',neuronSRMt3.triframes{1});
h5write('msID_SRMt3.h5','/frames/Test3littermates',neuronSRMt3.triframes{2});
h5write('msID_SRMt3.h5','/frames/Test3littermates_posnswap',neuronSRMt3.triframes{3});
h5write('msID_SRMt1.h5','/neuron/t',neuronSRMt1.t);
h5write('msID_SRMt2.h5','/neuron/t',neuronSRMt2.t);
h5write('msID_SRMt3.h5','/neuron/t',neuronSRMt3.t);


h5write('msID_SRMt1.h5','/neuron/IDarray',neuronSRMt1.IDarray);
h5write('msID_SRMt2.h5','/neuron/IDarray',neuronSRMt2.IDarray);
h5write('msID_SRMt3.h5','/neuron/IDarray',neuronSRMt3.IDarray);

