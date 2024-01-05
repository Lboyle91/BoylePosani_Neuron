function neuron=GetCNMFEData(filename);

neuron.name=filename;

neuron.YrA=h5read(filename,'/estimates/YrA');
neuron.C=h5read(filename,'/estimates/C');
neuron.S=h5read(filename,'/estimates/S');
neuron.ind=h5read(filename,'/estimates/A/indices');
neuron.indptr=h5read(filename,'/estimates/A/indptr');
neuron.A=h5read(filename,'/estimates/A/data');
neuron.dims=h5read(filename,'/estimates/dims');
neuron.idx_comp=h5read(filename,'/estimates/idx_components');
neuron.idx_comp_bad=h5read(filename,'/estimates/idx_components_bad');
neuron.bl=h5read(filename,'/estimates/bl');

sz=size(neuron.dims);

if sz(2)>2
    neuron.dims=[neuron.dims(1,1) neuron.dims(2,1)];
end

neuron.all_A(neuron.ind+1)=neuron.A;
neuron.nw_A=zeros(neuron.dims(1)*neuron.dims(2),length(neuron.indptr)-1);
neuron.all_A=zeros(neuron.dims(1)*neuron.dims(2),1);


for i=1:length(neuron.indptr)-1
    start=neuron.indptr(i)+1;
    last=neuron.indptr(i+1);
    neuron.nw_A(neuron.ind(start:last)+1,i)=neuron.A(start:last);
    
end

neuron.coor=zeros(length(neuron.indptr)-1,2);
for n=1:length(neuron.indptr)-1
    tempA=reshape(neuron.nw_A(:,n),neuron.dims(1),neuron.dims(2));
    high=0;
    for y=1:length(tempA(:,1))
       [pks,locs]=findpeaks(tempA(y,:));
                for p=1:length(pks)
                    if ~isempty(pks)
                        if pks(p)>high
                            high=pks(p);
                            neuron.coor(n,1)=locs(p);
                            neuron.coor(n,2)=y;
                        end
                    end
                end
    end
end
                           