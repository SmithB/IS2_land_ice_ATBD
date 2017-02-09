function write_ATL11_h5(h5_file, D3b, D3a);

if exist(h5_file,'file'); delete(h5_file); end


% write out the ATL06 data that comes with this file
ff=fieldnames(D3a);
for kf=1:length(ff);
    temp=D3a.(ff{kf});
    this_field=['/ATL06/', ff{kf}];
    h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
    h5write(h5_file, this_field,  temp);
end

% loop over the D3b structure
f0=fieldnames(D3b);
for k0=1:length(f0);
    temp0=D3b.(f0{k0});
    ff=fieldnames(temp0);
    % for each field in the D3b structure, create a field in a group
    for kf=1:length(ff);
        temp=temp0.(ff{kf});
        this_field=['/',f0{k0},'/',ff{kf}]; 
        if ~isempty(temp)
            h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
            h5write(h5_file, this_field,  temp);
        end
    end
end


% rewrite the directory structure
if false
    thedir='v10_hdf_subset/PIG_ATL11_6segs_dh_clouds_mat';
    [~, files]=unix(['ls ', thedir,'*/*/Track*.mat']);
    files=strsplit(deblank(files));
    for k=1:length(files);
        [thedir, thefile, extension]=fileparts(files{k}); 
        out_dir=strrep(thedir,'_mat','');
        if ~exist(out_dir,'dir'); mkdir(out_dir);end
        L=load(files{k});
        write_ATL11_h5([out_dir,'/', thefile,'.h5'], L.D_ATL11, L.D3);
    end


end
