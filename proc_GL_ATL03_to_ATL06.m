function proc_GL_ATL03_to_ATL06(sig_file, pair, this_rep)
 
% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})=h5read('SNR_F_table.h5', ['/',fields{kf}]);
end


 

D2_file=strrep(sig_file,'V052015.0_sigparms.h5','');
[D2a,  params]=read_ATLAS_GL_sim(D2_file, 1);
temp=regexp(  D2_file,'Track_(\d+)_','tokens');
this_RGT=str2double(temp{1}{1});

for k=1:2
    params(k).cycle=this_rep;
    params(k).RGT=this_RGT;
    params(k).PT=pair;
    params(k).GT=2*(params(k).PT-1)+k;
    params(k).orbit_number=(params(k).cycle-1)*1387+params(k).RGT;
    params(k).ATL03_sig_find=true;
end

[D3, dh_hist]=ATLAS_L3a_proc_ATBD(D2a,  params, [], SNR_F_table);
D3_file=strrep(D2_file,'.h5','_ATL06.mat');
save(D3_file,'D3','dh_hist','params');

if false
    sig_files=glob('/Volumes/ice1/ben/sdt/AnitaGroundFinder/GreenlandSimGF/TrackData_05/Pair_1/*parms.h5');
    fid=fopen('queue_GL_ATL03.txt','w');
    for k=1:length(sig_files)
        fprintf(fid, 'proc_GL_ATL03_to_ATL06(''%s'', 1, 5);\n', sig_files{k});
    end
    fclose(fid)
end


