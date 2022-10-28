clear D D_hist sigma_hat

% pick a big-but-reasonable window
Htot=30;
% pick Ntot based on a BGR of 8 MHz 
Ntot=ceil(57*Htot/1.5e8*8e6);
% Nsig is equivalent to between ~0.1 and ~2 PE/shot

sigma_vals=[0.1 0.25 0.5 1];
Nsig=0:10:Ntot-10;
reps=16;

out_dir='/home/ben/Dropbox/projects/IS2_ATBD/robust_peak_width_test_for_GSFC_July_2018';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

out_file=sprintf('%s/Ntot=%d_Htot=%d.h5', out_dir, Ntot, Htot);

for kN=1:length(Nsig)
    for ksig=1:length(sigma_vals)
        S=Nsig(kN);
        N=Ntot-S;
        for ii=1:reps
            this_S=max(0, S+round(randn(1)*sqrt(S)));
            this_N=max(0, N+round(randn(1)*sqrt(N)));
            D=[sigma_vals(ksig)*randn(this_S,1); Htot*(rand(this_N, 1)-0.5)];
            D=D(abs(D)<Htot/2);
            D_hist{ii,ksig, kN}=D;
            sigma_hat(ii, ksig, kN) =robust_peak_width_CDF(D, length(D)-S, [-0.5 0.5]*Htot);
            sigma_all(ii, ksig, kN)=iqr(D)/2;
        end
        med(ksig, kN)=median(sigma_hat(:, ksig, kN));
        med_all(ksig, kN)=median(sigma_all(:, ksig, kN));
        spread(ksig, kN)=iqr(sigma_hat(:, ksig, kN))/2;
    end
end

if true  
    outfile=sprintf('%s//Ntot=%d_Htot=%d.h5',out_dir, Ntot, Htot);
    if exist(outfile,'file'); delete(outfile); end
    for kN=1:length(Nsig)
        group1=sprintf('/Nsig_%d', Nsig(kN));
        for ksig=1:length(sigma_vals)
            group2=sprintf('%s/sigma_%2.2f', group1, sigma_vals(ksig));
            for ii=1:reps(end)
                group3=sprintf('%s/set%d', group2, ii)
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
    
    h5_file='~/Dropbox/projects/IS2_ATBD/robust_peak_width_test/Ntot=92_Htot=30.h5';
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
