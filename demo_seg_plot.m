function D3a=demo_seg_plot(ATL03_file, x_atc_0)

%files=glob('/Volumes/ice1/ben/sdt/ATLxx_example/PIG_Collab_v13B_NoFirn_NoDz_noDx/ATL03_subset/run_1/rep_*/Track_530-Pair_1_D2a.h5');

% note: for Track 530 Pair 1, rep 6 is great, rep 7 is so-so, 12 is pretty bad,
% smooth region;  2.8473e+07, +- 1 km
% rough region:  2.8470e+07 +- 1.5 km

[D2a, PairData, params, TrackData]=read_ATLAS_h5_D2a(ATL03_file, true);

if ~isfield(D2a,'seg_num')
    for kB=1:2
        D2a(kB).seg_num=1+floor(D2a(kB).x_RGT/20);
    end
end



ATL06_file=strrep(ATL03_file,'ATL03','ATL06');
ATL06_file=strrep(ATL06_file, 'D2.h5','D3.h5');
%[D3, dh_hist]=read_ATL06_h5(ATL06_file);


% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})=h5read('SNR_F_table.h5', ['/',fields{kf}]);
end

[~, this_seg_ind]=min(abs(D2a(2).x_RGT-x_atc_0));
this_seg=D2a(2).seg_num(this_seg_ind);

[D3a, dh_hist]=ATLAS_L3a_proc_ATBD('ATL03_to_ATL06',D2a, params, this_seg, SNR_F_table)
 