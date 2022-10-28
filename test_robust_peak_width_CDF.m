clear D D_hist sigma_hat

% pick a big-but-reasonable window
Htot=40;
% pick Ntot based on a BGR of 8 MHz 
Ntot=ceil(57*Htot/1.5e8*8e6/5)*5;  %N.B.  this is rounded up to the next value of 5
% pick Ntot based on a BGR of 12 MHz 
%Ntot=ceil(57*Htot/1.5e8*12e6);
REGEN=true
% Nsig is equivalent to between ~0.1 and ~2 PE/shot

sigma_vals=[0.1 0.25 0.5 1];
Nsig=0:5:Ntot;
reps=10240;

out_dir='/home/ben/Dropbox/projects/IS2_ATBD/robust_peak_width_test';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

if REGEN
    out_file=sprintf('%s/Ntot=%d_Htot=%d.h5', out_dir, Ntot, Htot);
    clear med med_all spread med_sigma
    for kN=1:length(Nsig)
        for ksig=1:length(sigma_vals)
            S=Nsig(kN);
            N=Ntot-S;
            %N=Ntot;%;-S;
            for ii=1:reps
                this_S=max(0, S+round(randn(1)*sqrt(S)));
                this_N=max(0, N+round(randn(1)*sqrt(N)));
                D=[sigma_vals(ksig)*randn(this_S,1); Htot*(rand(this_N, 1)-0.5)];
                D=D(abs(D)<Htot/2);
                D_hist{ii,ksig, kN}=D;
                sigma_hat(ii, ksig, kN) =robust_peak_width_CDF(D, length(D)-S, [-0.5 0.5]*Htot);
                sigma_all(ii, ksig, kN)=iqr(D, 0.5)/.6745;
                sigma(ii, ksig, kN)=std(D);
            end
            med_sigma(ksig, kN)=median(sigma(:, ksig, kN));
            med(ksig, kN)=median(sigma_hat(:, ksig, kN));
            med_all(ksig, kN)=median(sigma_all(:, ksig, kN));
            spread(ksig, kN)=iqr(sigma_hat(:, ksig, kN), 0.5)/.6745;
        end
    end
    save ~/Dropbox/projects/IS2_ATBD/robust_peak_width_test/Ntot=122_Htot=40_REVISED_stats Nsig sigma_vals med med_all med_sigma spread
end
if false 
    outfile=sprintf('%s//Ntot=%d_Htot=%d_REVISED.h5',out_dir, Ntot, Htot);
    if exist(outfile,'file'); delete(outfile); end
    for kN=1:length(Nsig)
        group1=sprintf('/Nsig_%d', Nsig(kN));
        for ksig=1:length(sigma_vals)
            group2=sprintf('%s/sigma_%2.2f', group1, sigma_vals(ksig));
            for ii=1:reps(end)
                group3=sprintf('%s/set%d', group2, ii);
                h5create(outfile,[group3,'/z'], [length(D_hist{ii, ksig, kN}), 1], 'Datatype','double');
                h5write(outfile,[group3,'/z'],D_hist{ii, ksig, kN});
                h5create(outfile,[group3,'/sigma_hat'], [1], 'Datatype','double');
                h5write(outfile,[group3,'/sigma_hat'],sigma_hat(ii, ksig, kN));
            end
        end
    end  
end

if false
    Htot=30;
    Ntot=ceil(57*Htot/1.5e8*8e6);
    sigma_vals=[0.1 0.25 0.5 1];
    Nsig=0:10:Ntot-10;
    reps=1:4;
    
    h5_file='~/Dropbox/projects/IS2_ATBD/robust_peak_width_test/Ntot=92_Htot=30_REVISED.h5';
    for kN=1:length(Nsig)
        group1=sprintf('/Nsig_%d', Nsig(kN));
        for ksig=1:length(sigma_vals)
            S=Nsig(kN);
            N=Ntot-S;
            group2=sprintf('%s/sigma_%2.2f', group1, sigma_vals(ksig));
            for ii=1:length(reps)
                group3=sprintf('%s/set%d', group2, ii);
                z=h5read(h5_file,[group3,'/z']);
                D_hist{ii,ksig, kN}=z;

                sigma_hat(ii, ksig, kN) =robust_peak_width_CDF(z, length(z)-S, [-0.5 0.5]*Htot);
                sigma_file(ii, ksig, kN)=h5read(h5_file,[group3,'/sigma_hat']);
            end
        end
    end
end

% Plot!
figure; clf; set(gcf,'units','inches','position', [1 1 8 4],'inverthardcopy','off','color','w','papersize', [8.5 4.5],'defaultaxesfontsize', 12);
hax=cheek_by_jowl(1, 4, [0.1 0.2 0.8 0.7]);
axes(hax(1)); imagesc(Nsig, 1:length(sigma_vals), med_sigma); caxis([0 24]); h_cb(1)=colorbar('northoutside');  xlabel(h_cb(1),'median(STD(z))')
axes(hax(2)); imagesc(Nsig, 1:length(sigma_vals), med_all); caxis([0 24]); h_cb(2)=colorbar('northoutside');  xlabel(h_cb(2),'median(\sigma_r(z))')
axes(hax(3)); imagesc(Nsig, 1:length(sigma_vals), med);  caxis([0 1.25]);  h_cb(3)=colorbar('northoutside');  xlabel(h_cb(3),'median(\sigma_{r,corr}(z))')
axes(hax(4)); imagesc(Nsig, 1:length(sigma_vals), spread); caxis([0 3.5]); h_cb(4)=colorbar('northoutside');  xlabel(h_cb(4),'spread(\sigma_{r,corr}(z))')

set(hax,'ytick', 1:length(sigma_vals),'yticklabel', num2str(sigma_vals(:)),'ydir','normal');

labels={'a','b','c','d'};
for k=1:4
    set(hax(k),'yticklabel', num2str(sigma_vals(:))); 
     ht(k)=text(0, 4.5, labels{k},'parent', hax(k));
     xlabel(hax(k),'N_{signal}');
end
set(ht,'horizontalalignment','left','verticalalignment','top','fontsize', 14,'backgroundcolor','w')

for k=2:3
    set(hax(k),'yticklabel','')
end
set(hax(4),'yaxislocation','right');


for k=[1 4]
    ylabel(hax(k),'\sigma_{signal}');
end

colormap(my_rgb_cpt(128).*repmat(linspace(0.5, 1, 128)', [1, 3]))
 
%colormap(my_rgb_cpt(128).*repmat([linspace(0.25, 0.75, 32)'; linspace(0.75, 1, 64+32)'], [1, 3]))




