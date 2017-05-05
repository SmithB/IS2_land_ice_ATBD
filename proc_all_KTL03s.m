

thedir='/Volumes/ice1/ben/sdt/KTL03/';
files=dir([thedir,'/ATL*.h5']);
out_dir=[thedir,'/ATL06'];
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

for kf=1:length(files)
    for kP=1:3
        out_file=sprintf('%s/ATL06/%s_Pair_%d.mat', thedir, strrep(files(kf).name,'.h5',''), kP);
        disp(out_file);
        beams=[2*kP-1 2*kP];
        try
            [H, ~, params]=read_sim_ATL03([thedir,'/', files(kf).name], kP);
        catch
            continue
        end
        PCE_vals=unique(cat(1, H.pce_mframe_cnt)); 
        PCE_vals=unique([1:100:max(PCE_vals), max(PCE_vals)]);
        
        clear D3a;
        all_segment_IDs=cat(1, H.ph_seg_num);
        all_PCE_vals=cat(1, H.pce_mframe_cnt);
        for kk=1:length(PCE_vals)-1
 
            for kB=1:2
                H1(kB)=index_struct(H(beams(kB)), ismember(H(beams(kB)).pce_mframe_cnt, max(1,PCE_vals(kk)-1):min(max(PCE_vals), PCE_vals(kk+1)+1)));
            end
            seg_id_list=unique(all_segment_IDs(ismember(all_PCE_vals, PCE_vals(kk):PCE_vals(kk+1)-1)));
            D3a1=KTL03_to_ATL06(H1, params(beams));
            D3a1=index_struct(D3a1, any(ismember(D3a1.segment_id, seg_id_list), 2)); 
            D3a(kk)=D3a1;
        end
        D3a=flatten_struct(D3a);
        save(out_file, 'D3a');
    end
end



% get the plots up on screen
ATL06_dir='/Volumes/ice1/ben/sdt/KTL03/ATL06/';
ATL03_dir='/Volumes/ice1/ben/sdt/KTL03/';
out_files=dir([ATL06_dir,'/*.mat']);
for kf=1:length(out_files)
    load([ATL06_dir,out_files(kf).name]);
    D6=D3a;
         
    YR=round_to(range(D6.h_med(isfinite(D6.h_med(:,1)) & D6.SNR_significance(:,1)>0.95)), 1)+[-10 10];
    
    
    S=regexp(out_files(kf).name,'(?<fname>\S+)_Pair_(?<pair>\d+).mat','names');
    D3=read_sim_ATL03([ATL03_dir, S.fname,'.h5'], str2double(S.pair));
    
    figure(kf);
    for ii=1:2
        kB=(str2num(S.pair)-1)*2+ii;
        subplot(2,1,ii);
        [II.z, II.x, II.y]=point_count_image(D3(kB).x_RGT, double(D3(kB).h_ph), 10, [YR(1) 1 YR(2)]);
        imagesc(II.x, II.y, II.z); hold on; caxis([0 100]); colormap(gray);
        plot(D6.x_RGT(:,kB), D6.h_mean(:, kB),'r.');
        els=D6.SNR_significance(:, kB)>0.99; 
        plot(D6.x_RGT(els,kB), D6.h_mean(els, kB),'bo');
      
    end
end






