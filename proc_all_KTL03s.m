

thedir='/Volumes/ice1/ben/sdt/KTL03/';
files=dir([thedir,'/ATL*.h5']);
out_dir=[thedir,'/ATL06'];
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})=h5read('SNR_F_table.h5', ['/',fields{kf}]);
end

for kf=1:length(files)
    for kP=1:3
        out_file=sprintf('%s/ATL06/%s_Pair_%d.mat', thedir, strrep(files(kf).name,'.h5',''), kP);
        disp(out_file);
        beams=[2*kP-1 2*kP];
        try
            [H, ~, params, dist_for_segment]=read_sim_ATL03([thedir,'/', files(kf).name], kP);
        catch
            continue
        end
        
        [params.orbit_number]=deal(1);
        [params.sigma_x]=deal(5);
        [params.t_dead]=deal(3.3e-9);
        for kB=1:2
            H(beams(kB)).track=ones(size(H(beams(kB)).x_RGT));
        end
        
        PCE_vals=unique(cat(1, H.pce_mframe_cnt)); 
        PCE_vals=unique([1:100:max(PCE_vals), max(PCE_vals)]);
        
        clear D3a;
        all_segment_IDs=cat(1, H.ph_seg_num);
        all_PCE_vals=cat(1, H.pce_mframe_cnt);
        for kk=1:length(PCE_vals)-1
            clear H1;
            for kB=1:2
                H1(kB)=index_struct(H(beams(kB)), ismember(H(beams(kB)).pce_mframe_cnt, max(1,PCE_vals(kk)-1):min(max(PCE_vals), PCE_vals(kk+1)+1)));
            end
            seg_id_list=unique(all_segment_IDs(ismember(all_PCE_vals, PCE_vals(kk):PCE_vals(kk+1)-1)));
            [D3a1, dh_hist1]=ATLAS_L3a_proc_ATBD(H1, params(beams), [], SNR_F_table, dist_for_segment(beams));
            %D3a1=KTL03_to_ATL06(H1, params(beams));
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
out_files=strsplit(deblank(ls([ATL06_dir,'/*version4*.mat'])));

colors=lines(2);

for kF=1:length(out_files)
    load(out_files{kF}); 
    D6=D3a;
    figure(kF); clf; 
    
    % Scale by 1000;
    D3a.x_RGT=(D3a.x_RGT-1.2441e7)/1000;
    D3a.dh_fit_dx=D3a.dh_fit_dx*1000;
    W_seg=40/1000;
   
    hax(1)=axes('position', [0.1 0.41, 0.8, 0.5]);
    hax(2)=axes('position', [0.1 0.1 0.8 0.3]);
    hold on; 
    for kB=1:2
        axes(hax(1)); hold on;
        hp=patch((D3a.x_RGT(:,kB)+[-1 1 1 -1]*W_seg/2)', (repmat(D3a.h_med(:,kB), [1, 4])+[0.5 0.5 -0.5 -0.5].*D3a.w_surface_window_final(:,kB) +W_seg/2*[-1 1 1 -1].*D3a.dh_fit_dx(:,kB))', colors(kB,:),'facealpha', 0.5);       
        set(hp,'edgecolor','none','facealpha', 0.5);
        good=D3a.SNR_significance(:, kB)<0.01;
        hl=plot_segs(D3a.x_RGT(good, kB), D3a.h_med(good, kB), D3a.dh_fit_dx(good, kB), W_seg,'k-');
        set(hl(ishandle(hl)),'linewidth', 2, 'color', colors(kB,:));
        plot(D3a.x_RGT(good, kB), D3a.h_med(good, kB),'o','color', colors(kB,:))
        plot(D3a.x_RGT(good, kB), D3a.h_med(good, kB),'k.','markersize',6);
        
        axes(hax(2)); hold on;
        plot(D3a.x_RGT(:, kB), D3a.SNR(:, kB),'.','color', colors(kB,:));
        plot(D3a.x_RGT(good, kB), D3a.SNR(good, kB),'o','color', colors(kB,:));
        plot(D3a.x_RGT(good, kB), D3a.SNR(good, kB),'k.', 'markersize', 6);
    end
    set(hax(1),'ylim', [-200 1e4]);
    set(hax(2),'yscale','log');
    linkaxes(hax,'x');
    set(hax(1),'xticklabel',[],'xgrid','on','ygrid','on')
    set(hax(2),'xgrid','on','ygrid','on')
    xlabel(hax(2),'along-track distance, m');
    ylabel(hax(1),'elevation'); ylabel(hax(2),'SNR');

end



% the current problem:
if false;
    XR=[1.2668    1.2673]*1e7;
    [H, ~, params, dist_for_segment]=read_sim_ATL03( '/Volumes/ice1/ben/sdt/KTL03//ATL03_for_ATL06_version4_sim1f_933_01.h5', kP);
    H2p=index_struct(H(2), H(2).x_RGT > XR(1) & H(2).x_RGT<XR(2));
    noise=H2p.h_ph < -45 | H2p.h_ph > -35;
    BGR=sum(noise)/(20*length(unique(double(H2p.ph_id_pulse)+double(H2p.pce_mframe_cnt*200)))/1.5e8)
    range(H2p.BGR)
end

if false
    for kf=1:length(out_files)
        load([ATL06_dir,out_files(kf).name]);
        D6=D3a;
        
        YR=round_to(range(D6.h_med(isfinite(D6.h_med(:,1)) & D6.SNR_significance(:,1)<0.01)), 1)+[-10 10];
        
        
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
    
end




