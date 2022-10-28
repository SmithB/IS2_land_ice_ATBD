WF.t=-3e-9:(2.5e-11):3e-9;
WF.p=gaussian(WF.t, 0, 0.68e-9);
 

thedir='~/Dropbox/scf/';
files=dir([thedir,'/ATL03*.h5']);
out_dir=[thedir,'/ATL06_matlab'];
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% read in the SNR F table:
fields={'BGR', 'W_surface_window_initial','SNR', 'P_NoiseOnly'};
for kf=1:length(fields)
    SNR_F_table.(fields{kf})=h5read('SNR_F_table_june26_2018.h5', ['/',fields{kf}]);
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
        %[params.WF]=deal(WF);
        for kB=1:2
            H(beams(kB)).track=ones(size(H(beams(kB)).x_RGT));
        end
        
        PCE_vals=unique(cat(1, H.pce_mframe_cnt)); 
        PCE_vals=unique([1:100:max(PCE_vals), max(PCE_vals)]);
        
        clear D3a;
        [D3a, dh_hist]=ATLAS_L3a_proc_ATBD(H(1:2), params(beams), [], SNR_F_table, dist_for_segment(beams));
        figure; h=cheek_by_jowl(3, 1, [0.05 0.05 0.9 0.9]);
        for kB=1:2; 
            axes(h(kB)); 
            plot(H(kB).x_RGT, H(kB).h_ph,'.','markersize', 1); hold on;
            plot_segs(D3a.x_RGT(:,1), D3a.h_mean(:,1), D3a.dh_fit_dx(:,1),40,'r');
        end
        axes(h(3));hold on;
        for kB=1:2
            plot(D3a.x_RGT(:,kB), D3a.n_fit_photons(:, kB),'o');
        end
                
        save(out_file, 'D3a');
    end
end



% get the plots up on screen
ATL03_dir='/Volumes/ice1/ben/sdt/Alaska_sim/';
ATL06_dir=[ATL03_dir,'/ATL06'];
 
out_files=strsplit(deblank(ls([ATL06_dir,'/*Pair_1.mat'])));

colors=lines(2);

for kF=1:length(out_files)
    load(out_files{kF}); 
    D6=D3a;
    figure(kF); clf; 
    
    % Scale by 1000;
    D3a.x_RGT=(D3a.x_RGT-13895760)/1000;
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
    set(hax(1),'ylim', [200 1400]);
    set(hax(2),'yscale','log');
    linkaxes(hax,'x');
    set(hax(1),'xticklabel',[],'xgrid','on','ygrid','on')
    set(hax(2),'xgrid','on','ygrid','on')
    xlabel(hax(2),'along-track distance, m');
    ylabel(hax(1),'elevation'); ylabel(hax(2),'SNR');
    title(hax(1),strrep(out_files(kF),'_','-'));
end

 

