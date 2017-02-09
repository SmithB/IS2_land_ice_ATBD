function D=read_matlas(filename)


pair={'left','center','right'};

fields={'delta_time','latitude','longitude','height', 'photon_type','shot_num'};
for kP=1:length(pair)
    GG=['/photon/', pair{kP}];
    channels=h5info(filename,GG);
    channels={channels.Groups.Name};
    for kC=1:length(channels)
        for kf=1:length(fields)
            temp.(fields{kf})=h5read(filename, [channels{kC},'/strong/', fields{kf}]); 
        end
        temp.x_RGT=temp.delta_time*200;
        temp.time=temp.delta_time/24/3600;
        temp.h=temp.height;
        temp=rmfield(rmfield(temp, 'delta_time'),'height');
        rr=regexp(channels{kC},'channel(\d+)','tokens');
        temp.channel=zeros(size(temp.time))+str2double(rr{1}{1});
        D.(pair{kP})(kC)=temp;
    end

end
    
    
   

 