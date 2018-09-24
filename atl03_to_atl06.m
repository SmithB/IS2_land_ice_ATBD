function [D_06, dh_hist]=atl03_to_atl06(fname, pair, SNR_F_table_file)
 
 
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

[D_06, dh_hist]=deal(cell(length(pair),1));
for kP=pair
    beams=[2*kP-1 2*kP];
    [H, ~, params, dist_for_segment]=read_ATL03(fname, kP);
    [params.orbit_number]=deal(1);
    [params.sigma_x]=deal(5);
    [params.t_dead]=deal(3.3e-9);
    %[params.WF]=deal(WF);
    for kB=1:2
        H(beams(kB)).track=ones(size(H(beams(kB)).x_RGT));
    end
    
    PCE_vals=unique(cat(1, H.pce_mframe_cnt));
    PCE_vals=unique([1:100:max(PCE_vals), max(PCE_vals)]);
    
    clear D3a;
    [D_06{kP}, dh_hist{kP}]=ATLAS_L3a_proc_ATBD(H(1:2), params(beams), [], SNR_F_table, dist_for_segment(beams));
end


