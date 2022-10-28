function [D_06, dh_hist, D_03]=atl03_to_atl06(fname, pair, SNR_F_table_file, out_dir, which_seg)
 
if ~exist('pair','var') || isempty(pair) 
    pair=[1,2,3];
end

if ~exist('SNR_F_table_file','file') || isempty(SNR_F_table_file)
    % assume that SNR_f_table lives in thesame directory as the
    % atl03_to_atl06_script
    [thedir, ~]=fileparts(which(mfilename));
    SNR_F_table_file= [thedir,'/','SNR_F_table_June26_2018.h5'];
end

% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})=h5read(SNR_F_table_file, ['/',fields{kf}]);
end

if ~exist('out_dir','var') || isempty(out_dir)
    out_dir=pwd;
end

% make sure the file is valid
beams={'l','r'};
valid_beams=true(2,3);
for kP=1:3
    for kB=1:2
        try
            II=h5info(fname, sprintf('/gt%d%s', kP, beams{kB}));
        catch
            valid_beams(kB, kP)=false;
        end
    end
end

if ~any(all(valid_beams))
    D_06=[]; dh_hist=[]; D_03=[];
    return
end

D_06=cell(length(pair),1);
Dh_hist=cell(2, length(pair));
for kP=pair(:)'
    if ~all(valid_beams(:, kP))
        continue
    end
    beams=[2*kP-1 2*kP];
    Geoloc_full=read_ATL03_geolocation(fname, kP, [], true);
    params=read_ATL03_params(fname, kP);
    [params.orbit_number]=deal(1);
    [params.sigma_x]=deal(5);
    [params.t_dead]=deal(3.2e-9);
    
    [uID]=unique(cat(1, Geoloc_full.segment_id));
    if length(uID) >=500
        ID0=uID(1):500:uID(end); ID0(end)=uID(end);
    else
        ID0=[uID(1), uID(end)];
    end
    
    if exist('which_seg','var') && ~isempty(which_seg);
        ID0=ID0(abs(ID0-which_seg)<500);
    end
    
    clear D6a dh_hist_a
    %disp('WARNING: ONLY CALCULATING 500 segments');
    for k0=1:length(ID0)-1
        this_seg_range=[max(ID0(k0)-10, 1), min(ID0(k0+1)+20, ID0(end))];
        proc_segs=ID0(k0):ID0(k0+1)+9;
        output_segs=ID0(k0):ID0(k0+1)-1;
        [geoloc_sub, dist_for_segment]=read_ATL03_geolocation(fname, kP,  struct('seg_range', this_seg_range), Geoloc_full);
        if isempty(geoloc_sub)|| (isempty(geoloc_sub(beams(1)).delta_time) & isempty(geoloc_sub(beams(2)).delta_time)); continue; end
        D_03=read_ATL03_photon_data(fname, kP, [], geoloc_sub);
        for kB=1:2
            D_03(kB).BGR(~isfinite(D_03(kB).BGR))=1e6;
        end
        if exist('which_seg','var') && ~isempty(which_seg);
            for kB=1:2
                D_03(kB)=index_struct(D_03(kB), D_03(kB).segment_id==which_seg | D_03(kB).segment_id==which_seg-1);
                proc_segs=proc_segs(proc_segs==which_seg);
            end
        end
        
        [temp1, temp2]=ATLAS_L3a_proc_ATBD(D_03(beams), params(beams), proc_segs, SNR_F_table, dist_for_segment(beams));
        if ~isempty(temp1)
            for kB=1:2
                % select the members of
                cols=ismember(min(temp2(kB).segment_id_list), output_segs);
                for ff=fieldnames(temp2(kB))'
                    if ~strcmp(ff{1},'dh')
                        temp2(kB).(ff{1})=temp2(kB).(ff{1})(:, cols);
                    end
                end
            end
            temp1=index_struct(temp1, any(ismember(temp1.seg_count, output_segs), 2));
            
            [D6a(k0), dh_hist_a(k0,:)]=deal(temp1, temp2);
        end
        if mod(k0, 20)==0;
            fprintf(1, 'k0=%f out of %f\n', k0, length(ID0));
        end
    end
    for field = fieldnames(D6a)'
        D_06{kP}.(field{1})=cat(1, D6a.(field{1}));
    end
    for field = fieldnames(dh_hist_a)'
        for kB=1:2
            dh_hist{kB, kP}.(field{1})=cat(2, dh_hist_a(:, kB).(field{1}));
        end
    end
end

[~, thefile, ~]=fileparts(fname);

versions={'_943_','_944_','_200_','_201_'};
for kv=1:length(versions)
    thefile=strrep(thefile, versions{kv},'_mat_');
end

out_file=[out_dir,'/',strrep(thefile,'ATL03','ATL06'),'.h5'];
ATL06_structure_gen(cat(1, D_06{:}), dh_hist, out_file);

    

