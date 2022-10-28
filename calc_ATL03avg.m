function calc_ATL03avg(ATL03_file, out_dir)

[~, file_base]=fileparts(ATL03_file);

out_file=[out_dir,'/',strrep(file_base,'ATL03_', 'ATL03Avg_'),'.h5'];
disp(out_file)
if exist(out_file,'file')
    return
end
for pair=1:3
    try
        D03=read_ATL03_photon_data(ATL03_file, pair);
    catch
        continue
    end
    fprintf(1,'working on %s pair %d\n', ATL03_file, pair);
    D03=D03(2*(pair-1)+[1 2]);
    for kB=1:2
        D03(kB)=index_struct(D03(kB), D03(kB).lat_ph < -87.5 & D03(kB).signal_conf_ph >= 1);
    end
    uSeg=unique([D03(1).segment_id; D03(2).segment_id]);
    [S.hbar_ge1, S.hbar_ge3_plus, S.hbar_ge3,  S.sigma_ge1, S.sigma_ge3, S.hbar_TSE, S.sigma_TSE, latbar, lonbar]=deal(NaN(length(uSeg),2));
    
    for kB=1:2
        elsB=bin_by(D03(kB).segment_id, uSeg);
        for k_bin=2:length(elsB)
            bins=(uSeg==uSeg(k_bin) | uSeg==uSeg(k_bin)-1);
            if ~any(bins); continue; end
            ii=cat(1, elsB{bins});
            if length(ii) > 20
                ss=D03(kB).signal_conf_ph(ii);
                hh=D03(kB).h_ph(ii);
                S.hbar_ge1(k_bin, kB)=mean(hh);
                S.sigma_ge1(k_bin, kB)=std(D03(kB).h_ph(ii));

                ii_sub=ss >=3;
                if sum(ii_sub) > 20
                    S.hbar_ge3(k_bin, kB)=mean(hh(ii_sub));
                    S.sigma_ge3(k_bin, kB)=std(hh(ii_sub));
                
                    S.hbar_ge3_plus(k_bin, kB)=mean(hh(abs(hh-S.hbar_ge3(k_bin, kB))<1.5));
                    [S.hbar_TSE(k_bin, kB), S.sigma_TSE(k_bin, kB)]=TSE_h(D03(kB).x_RGT(ii), hh, ss);
                end
                latbar(k_bin, kB)=mean(D03(kB).lat_ph(ii));
                lonbar(k_bin, kB)=mean(D03(kB).lon_ph(ii));
            end
        end
    end
    if length(uSeg)==0
        continue
    end
    DS=sprintf('/GT%d/segID', pair);
    h5create(out_file,DS, length(uSeg));
    h5write(out_file, DS, uSeg)
    
    for field={ 'hbar_ge1', 'hbar_ge3_plus', 'hbar_ge3',  'sigma_ge1', 'sigma_ge3', 'hbar_TSE', 'sigma_TSE'}
        
        DS=sprintf('/GT%d/%s', pair, field{1});
        h5create(out_file, DS, [length(uSeg), 2]);
        h5write(out_file, DS, S.(field{1}));
    end
   
    DS=sprintf('/GT%d/latitude', pair);
    h5create(out_file, DS, [length(uSeg), 2]);
    h5write(out_file, DS, latbar);
    
    DS=sprintf('/GT%d/longitude', pair);
    h5create(out_file, DS, [length(uSeg), 2]);
    h5write(out_file, DS, lonbar);
    
end

%--------------------------------
function [h0, sigma]=TSE_h(x, h, s)

x=x-mean(range(x));
G=[ones(size(x)), x];
good=s>=3;
last_good=false(length(x),1);
HW=10;
h0=mean(h(good));
count=0;
while count<10 && sum(abs(last_good-good))>1 
    count=count+1;
    last_good=good;
    G1=G(good,:);
    m=G1\h(good);
    r=h-G*m;
    sigma=iqr(r(good))/2;
    HW=max([3, 6*sigma, 0.75*HW]);
    good=abs(r)<HW/2;
end
h0=m(1);



if false
    fid=fopen('calc_ATL03_avg_queue.txt','w');
    files=glob('ATL03*.h5');
    for k=1:length(files);
        fprintf(fid,'calc_ATL03avg(''%s'',''ATL03_avg'');\n', files{k});
    end
end
       



