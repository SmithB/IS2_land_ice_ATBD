function dump_struct_to_h5(IN, out_file, groupname, template)  

if ~exist('template','var')
    template={'E%d'};
end

if length(template)==1
    template={template{1},'E%d'};
end
if isstruct(IN) & length(IN(:))>1
    for k=1:length(IN(:))
        if min(size(IN))>1
            [r,c]=ind2sub(size(IN), k);
            dump_struct_to_h5(IN(k), out_file, [groupname,'/', sprintf(template{1}, r, c)],template(2:end));
        else
            dump_struct_to_h5(IN(k), out_file, [groupname,'/', sprintf(template{1}, k)],template(2:end));
        end
    end
elseif isstruct(IN)
    fields=fieldnames(IN);
    for kf=1:length(fields)
        ff=fields{kf};
        fc=IN.(ff);
        if isempty(fc); continue; end
        if ~isstruct(fc)
            DSname=[groupname,'/', ff];
            if islogical(fc)
                h5create(out_file, DSname, size(fc),'datatype','int8','chunksize', size(fc));
                h5write(out_file, DSname, int8(fc));
            else
                h5create(out_file, DSname, size(fc),'datatype','double','chunksize', size(fc));
                h5write(out_file, DSname, double(fc));
            end
        else
            dump_struct_to_h5(fc, out_file,[groupname,'/',ff], template(2:end));   
        end
    end
end
     
     
        
        
    


if false
    [~, files]=unix('ls ~/Dropbox/projects/IS2_ATBD/Combined_corr_test_data_Jul_7_2017/D3/*.mat');
    files=strsplit(deblank(files));
    %files={'/home/ben/Dropbox/projects/IS2_ATBD/Combined_corr_test_data_Jul_7_2017/D3/Rough=9.60e-01_Rsurf=1.25e-01_BGR=1.00e+07_LOG.mat'};
    for k0=1:length(files)
        clear LOG
        out_file=strrep(files{k0},'.mat','.h5');
        if exist(out_file,'file'); delete(out_file); end
        load(files{k0});
        dump_struct_to_h5(LOG, out_file,'',{'seg%d_beam%d','it_%d'});
    end 
end


