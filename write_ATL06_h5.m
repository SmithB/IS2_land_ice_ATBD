function write_ATL06_h5(D3, dh_hist, h5_file);


if exist(h5_file,'file'); delete(h5_file); end

ff=fieldnames(D3);
for kf=1:length(ff);
    temp=D3.(ff{kf});
    this_field=['/', ff{kf}];
    h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
    h5write(h5_file, this_field,  temp);
end

ff=fieldnames(dh_hist);
ff=ff(~ismember(ff, 'count'));
beam_names={'L/','R/'};
for kB=1:2;
    for kf=1:length(ff);
        this_field=['/hist_200m/', beam_names{kB},ff{kf},'/'];
        temp=dh_hist(kB).(ff{kf});
        temp=temp(:);
        h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
        h5write(h5_file, this_field,  temp);
    end
    this_field=['/hist_200m/', beam_names{kB},'count/'];
    temp=dh_hist(kB).count;
    h5create(h5_file, this_field, size(temp),'ChunkSize', [size(temp,1), min(size(temp,2), 20)], 'Datatype','uint16','Deflate', 9);
    h5write(h5_file, this_field, temp);
end