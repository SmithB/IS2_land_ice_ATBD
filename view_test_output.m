%test_dir='/home/ben/Dropbox/projects/IS2_ATBD/Combined_corr_test_data_Jul_25_2017/';
test_dir='/home/ben/Dropbox/projects/IS2_ATBD/Combined_corr_test_data/Nov_6_2017_all_SNR/';
[~,files]=unix(sprintf('ls %s/D3/Rough*07.h5', test_dir));
files=strsplit(deblank(files));

clear D
for k=1:length(files); 
    temp=regexp(files{k},'Rough=(\S+)_Rsurf=(\S+)_BGR=(\S+).h5','tokens');
    D(k).Rough=str2num(temp{1}{1}); 
    D(k).Rsurf=str2num(temp{1}{2});
    D(k).BGR=str2num(temp{1}{3});
end

fields={'w_surface_window_final','h_LI','dh_fit_dx'};
clear D1;
Rough0=unique([D.Rough]);
Rsurf0=unique([D.Rsurf]);
[Rough, Rsurf]=meshgrid(Rough0, Rsurf0);
for kc=1:size(Rough, 2)
 
    figure(kc);clf; hax=cheek_by_jowl(5,3, [0.15 0.15 0.7 0.7]);
    title(hax(1,1), sprintf('Rough=%f', Rough(1,kc)));
     for kr=1:size(Rough,1)
         k=find([D.Rough]==Rough(kc) & [D.Rsurf]==Rsurf(kr));
        for kf=1:length(fields) 
            V=h5read(files{k},['/', fields{kf}]);
            for beam=1:2
                temp=V(isfinite(V(:,beam)), beam);
                D1(beam).(fields{kf}).med(kr, kc)=median(temp);
                D1(beam).(fields{kf}).bar(kr, kc)=mean(temp);
                D1(beam).(fields{kf}).sigma_hat(kr, kc)=iqr(temp)/2;
                D1(beam).(fields{kf}).P95(kr, kc)=percentile(temp, 0.95);
                D1(beam).(fields{kf}).P05(kr, kc)=percentile(temp, 0.05);
            end
        end
        D3=read_ATL06_h5(files{k});
        axes(hax(kr,1)); hold on;
        plot_segs(D3.x_RGT(:,2), D3.h_LI(:,2), D3.dh_fit_dx(:,2),40,'r');
        plot_segs(D3.x_RGT(:,1), D3.h_LI(:,1), D3.dh_fit_dx(:,1),40,'b');
        ylabel(sprintf('R-surf=%4.3f', Rsurf(kr))); set(gca,'ylim', [-10 20]);
        
        axes(hax(kr,2)); hold on; 
        plot(D3.x_RGT(:,2), D3.SNR_significance(:,2),'r.');
        plot(D3.x_RGT(:,1), D3.SNR_significance(:,1),'bo');
        set(gca,'ylim', [0 1]);
        
        axes(hax(kr,3)); hold on; 
        plot(D3.x_RGT(:,2),D3.signal_selection_source(:,2),'r.');
        plot(D3.x_RGT(:,1),D3.signal_selection_source(:,1),'bo');
        set(gca,   'ylim', [-0.5 2.5]);    
     end
 
end


