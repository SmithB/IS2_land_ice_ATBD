% ...quick test of BGR calculation...


ATL06_dir='/Volumes/ice1/ben/sdt/KTL03/ATL06/';
ATL03_dir='/Volumes/ice1/ben/sdt/KTL03/';
out_files=dir([ATL06_dir,'/*.mat']);
for kf=1:length(out_files)
    load([ATL06_dir,out_files(kf).name]);
    D6=D3a;     
    YR=round_to(range(D6.h_med(isfinite(D6.h_med(:,1)) & D6.SNR_significance(:,1)>0.95)), 1)+[-10 10];
    
    S=regexp(out_files(kf).name,'(?<fname>\S+)_Pair_(?<pair>\d+).mat','names');
    D3=read_sim_ATL03([ATL03_dir, S.fname,'.h5'], str2double(S.pair));
    
     
    figure(kf); clf;
    clear BG_check
    for kB=1:2
        [ee,uPCE]=bin_by(D3(kB).pce_mframe_cnt);
        uPCE=unique(D3(kB).pce_mframe_cnt);
        [dY, N, Nshots, BGR]=deal(NaN(size(uPCE)));
        for kP=1:length(uPCE)
            els=ee{kP};
            if length(els) > 10
                dY(kP)=diff(percentile(D3(kB).h_ph(els), [0.05 0.95]));
                N(kP)=length(els);
                Nshots(kP)=diff(range(D3(kB).pulse_num(els)));
                BGR(kP)=median(D3(kB).BGR(els));
            end
        end
        BG_check(kB)=struct('dY', dY,'N', N,'Nshots', Nshots,'BGR', BGR, 'N_obs',  N./(200.*dY/(1.5e8)), 'PCC', uPCE);
        subplot(2,2,kB);
        plot(uPCE, BG_check(kB).BGR/1e6,'b.'); hold on;
        plot(uPCE, BG_check(kB).N_obs/1e6,'r.'); 
        subplot(2,2, kB+2);
        histogram(BG_check(kB).N_obs./BG_check(kB).BGR, 0:.1:4)
    end
    save(['/Volumes/ice1/ben/sdt/KTL03/BGR_check/',out_files(kf).name],'BG_check') 
end
    




    