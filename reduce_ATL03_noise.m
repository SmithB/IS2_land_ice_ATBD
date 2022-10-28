function [D, params, TrackData, PairData]=reduce_ATL03_noise(in_h5,  target_BGR, out_h5  )

% [D, params, TrackData, PairData]=reduce_ATL03_noise(in_h5, target_BGR, out_h5)
% remove a fraction of the noise in ATL03 D2 file 'in_h5' to produce an
% equivalent data set, with noise rate equal to target_BGR.  
% target_BGR must be less that the background rate in in_h5.
% Optionall write out the dataset to file 'out_h5'

if iscell(in_h5)
     [D, PairData, params, TrackData]=deal(in_h5{1:4});
else
    [D, PairData, params, TrackData]=read_ATLAS_h5_D2a(in_h5, true);
end
for kB=1:length(D)
    ind=D(kB).SigNoise==1 |...
        ( D(kB).SigNoise==0 & rand(size(D(kB).SigNoise))<target_BGR./D(kB).BGR);
    D(kB)=index_struct(D(kB), ind);
    D(kB).BGR(:)=target_BGR;
    params(kB).NoiseRate=target_BGR;
end


if exist('out_h5','var')  
    [out_dir, thefile]=fileparts(out_h5);
    if ~exist(out_dir,'dir');
        mkdir(out_dir);
    end
    write_D2_HDF(struct('file', out_h5,'D2a', D, 'params',params, 'TrackData', TrackData,'PairData' , PairData), out_dir)
end

% run as a batch:
if false
    base='/Volumes/ice1/ben/sdt/ATLxx_example/v13_hdf/PIG_ATL03_8.00MHz';   
    BGR_out=[0.25 0.5 1 2 4]; % MHz;
    [~, out]=unix(sprintf('ls %s/tau=*/rep*/*D2.h5 | sort -r',base));
    ATL03_list=strsplit(deblank(out));
    for kF=1:length(ATL03_list)   
        if isempty(ATL03_list{kF}); continue; end
        out_files={};
        for kB=1:length(BGR_out)
            out_top=sprintf('PIG_ATL03_%3.2fMHz', BGR_out(kB));
            out_files{kB}=strrep(ATL03_list{kF},'PIG_ATL03_8.00MHz',out_top);
            DoThisOne(kB)=~exist(out_files{kB},'file');
        end
        
        out_files=out_files(DoThisOne);
        if isempty(out_files); continue; end
        clear temp;
        [temp{1}, temp{2}, temp{3}, temp{4}]=read_ATLAS_h5_D2a(ATL03_list{kF}, true);
        for kB= 1:length(BGR_out)
            if ~DoThisOne(kB); continue; end
            out_dir=fileparts(out_files{kB});
            if ~exist(out_dir,'dir'); mkdir(out_dir); end
            fprintf(1,'working on %s\n', out_files{kB});
            reduce_ATL03_noise( temp,  BGR_out(kB)*1e6, out_files{kB});
        end
    end
end