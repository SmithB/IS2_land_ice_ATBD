function [D_06, dh_hist, D_03]=atl03_to_atl06(fname, pair, SNR_F_table_file)
 
 
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


D_06=cell(length(pair),1);
Dh_hist=cell(2, length(pair));
for kP=pair(:)'
    beams=[2*kP-1 2*kP];
    Geoloc_full=read_ATL03_geolocation(fname, kP, [], true);
    params=read_ATL03_params(fname, beams);
    [params.orbit_number]=deal(1);
    [params.sigma_x]=deal(5);
    [params.t_dead]=deal(3.3e-9);
    
    [uID]=unique(cat(1, Geoloc_full.segment_id));
    ID0=uID(1):500:uID(end); ID0(end)=uID(end);
    
    clear D6a dh_hist_a
    for k0=1:length(ID0)-1
        [geoloc_sub, dist_for_segment]=read_ATL03_geolocation(fname, kP,  struct('seg_range', [ID0(k0), ID0(k0+1)-1]), Geoloc_full);
        if isempty(geoloc_sub)|| (isempty(geoloc_sub(beams(1)).delta_time) & isempty(geoloc_sub(beams(2)).delta_time)); continue; end
        D_03=read_ATL03_photon_data(fname, kP, [], geoloc_sub);
        %         for kB=1:2
        %             D_03(beams(kB)).track=ones(size(D_03(beams(kB)).x_RGT));
        %         end
        [temp1, temp2]=ATLAS_L3a_proc_ATBD(D_03(beams), params(beams), ID0(k0):ID0(k0+1)-1, SNR_F_table, dist_for_segment(beams));
        if ~isempty(temp1)
            [D6a(k0), dh_hist_a(k0,:)]=deal(temp1, temp2);
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


