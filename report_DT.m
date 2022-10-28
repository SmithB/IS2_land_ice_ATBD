 


data_dir_var='/home/ben/Dropbox/projects/IS2_ATBD/deadtime_correction_data/vdt/CAL_42_PRIM_2016298_a';
data_dir_const=[data_dir_var,'_meanDT'];

[~,files]=unix(['ls ', data_dir_var,'/D2']);
files=strsplit(deblank(files));

for k=1:length(files)
    temp=regexp(files{k},'Rough=(.*)_Rsu','tokens');
    Rough(k)=str2double(temp{1});
    D2=make_fpb_test_data_var_deadtime('read_test_file',  [data_dir_var,'/D2/', files{k}]);
    for kB=1:6    
        bar_v(k, kB)=mean(D2(kB).h(D2(kB).detected~=0));
        med_v(k, kB)=median(D2(kB).h(D2(kB).detected~=0));
        sigma_v(k, kB)=std(D2(kB).h(D2(kB).detected~=0))/sqrt(sum(D2(kB).detected~=0));
        %fprintf(1, 'beam=%d, mean=%f, median=%f\n', kB, mean(D2(kB).h(D2(kB).detected~=0)),median(D2(kB).h(D2(kB).detected~=0)));
    end
    
         
    D2=make_fpb_test_data_var_deadtime('read_test_file',  [data_dir_const,'/D2/', files{k}]);
    for kB=1:6    
        bar_c(k, kB)=mean(D2(kB).h(D2(kB).detected~=0));
        med_c(k, kB)=median(D2(kB).h(D2(kB).detected~=0));
        sigma_c(k, kB)=std(D2(kB).h(D2(kB).detected~=0))/sqrt(sum(D2(kB).detected~=0));
        %fprintf(1, 'beam=%d, mean=%f, median=%f\n', kB, mean(D2(kB).h(D2(kB).detected~=0)),median(D2(kB).h(D2(kB).detected~=0)));
    end
         
end

figure(1); clf; subplot(2,2,1); hold on; 
colors=lines(6);
[~, ind]=sort(Rough);
for kB=1:6
   HH(kB)=errorbar(Rough(ind), bar_v(ind, kB), sigma_v(ind, kB),  '-','color', colors(kB,:));
    errorbar(Rough(ind), bar_c(ind, kB), sigma_c(ind, kB),'--','color', colors(kB,:));
end
legend(HH, '1W','1S','2W','2S','3W','3S')
xlabel('roughness, m'); ylabel('mean bias, m');

subplot(2,2,3); hold on;
for kB=1:6
   errorbar(Rough(ind), bar_v(ind, kB)-bar_c(ind, kB), sqrt(sigma_v(ind, kB).^2+sigma_c(ind, kB).^2),'-','color', colors(kB,:));
end
xlabel('roughness, m'); ylabel('variable -constant mean bias, m');


subplot(2,2,2); hold on; 
colors=lines(6);
[~, ind]=sort(Rough);
for kB=1:6
    errorbar(Rough(ind), med_v(ind, kB),sigma_v(ind, kB),'-','color', colors(kB,:));
    errorbar(Rough(ind), med_c(ind, kB),sigma_c(ind, kB),'--','color', colors(kB,:));
end
xlabel('roughness, m'); ylabel('median bias, m');


subplot(2,2,4); hold on;
for kB=1:6
   errorbar(Rough(ind), med_v(ind, kB)-med_c(ind, kB), sqrt(sigma_v(ind, kB).^2+sigma_c(ind, kB).^2),'-','color', colors(kB,:));
end
xlabel('roughness, m'); ylabel('variable -constant median bias, m');
set(findobj('type','ErrorBar'),'linewidth', 2)


clear D;
D{1}=csvread('~/Dropbox/projects/IS2_ATBD/deadtime_by_channel/CAL_42_PRIM_2016298_a.csv', 17,0);
[~, out]=unix('head -11 ~/Dropbox/projects/IS2_ATBD/deadtime_by_channel/CAL_42_PRIM_2016298_a.csv | tail -1'); 
temp=regexp(out,',(.*),','tokens');
T(1)=str2num(temp{1}{1}); 

D{2}=csvread('~/Dropbox/projects/IS2_ATBD/deadtime_by_channel/CAL_42_PRIM_2016298_b.csv', 17,0);
[~, out]=unix('head -11 ~/Dropbox/projects/IS2_ATBD/deadtime_by_channel/CAL_42_PRIM_2016298_b.csv | tail -1'); 
temp=regexp(out,',(.*),','tokens');
T(2)=str2num(temp{1}{1}); 


D{3}=csvread('~/Dropbox/projects/IS2_ATBD/deadtime_by_channel/CAL_42_PRIM_2016298_c.csv', 17,0);
[~, out]=unix('head -11 ~/Dropbox/projects/IS2_ATBD/deadtime_by_channel/CAL_42_PRIM_2016298_c.csv | tail -1'); 
temp=regexp(out,',(.*),','tokens');
T(3)=str2num(temp{1}{1}); 
[~, iT]=sort(T); 

figure; clf;;
h=cheek_by_jowl(6, 10, [0 0 1 1])';
for k=1:60
    axes(h(k)); hold on;
    for kT=1:3
        DT(kT)= D{kT}(k, 3);
    end
    plot(T(iT), DT(iT),'marker','o')
end




