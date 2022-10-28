function [recorded, det, dead_all] = deadtime2(t, pulse,  N_det, t_dead_analog, t_dead_digital, gain_params)
% NOTE: t and pulse must be sorted by pulse, then by time.
 
det = ceil(rand(size(t))*N_det);
recorded = false(size(t));

if numel(t_dead_digital)==1
    % If there's just one value for t_dead_digital, replicate it over
    % channels and leading vs falling edge
    t_dead_digital=zeros(N_det, 2)+t_dead_digital;
end

t_dead_digital_max = max(t_dead_digital(:));
t_dead_digital_min = min(t_dead_digital(:));

% the analog detector limits the throughput to the digital detector so
% you only need to search a maximum of t_dead_digita/t_dead_analog
% photons away from a detection event

maxPhtn = ceil(t_dead_digital_max/t_dead_analog);
    
% loop for each detector
for i = 1:N_det
    
    idx = det == i;
    t0 = t(idx);
    pulse0 = pulse(idx);
    
    % ----------- Analog paralyzable detector -----------
    % photons that arive at least t_dead_analog after last event 
    dt0 = [0; (t0(2:end) -t0(1:(end-1)))];
    recorded0 = dt0 > t_dead_analog;
    
    % first photon in pulse always regesters
    firstPhtn = [true; (pulse0(2:end) -pulse0(1:(end-1))) ~= 0];
    
    % photon events that pass through to the Analog detector to the digital
    % detector
    idxA = recorded0 | firstPhtn;
    
    % ----------- Digital non-paralyzable electronics -----------
    
    % time of events that make it through to the digitalizing electronics
    t0 = t0(idxA);
    dt0 = dt0(idxA);
    firstPhtn = firstPhtn(idxA);
    
    % find recorded events assuming paralyzable electronics (first pass)
    recorded0 = (dt0 > t_dead_digital_max) | firstPhtn;
    maxSrch = length(t0);

    % only search photons event that are spaced < max(t_dead_digitalA,
    % t_dead_digitalB) to see if they fall after digital detector deadtime
    chkPhtn = ~recorded0;

    while any(chkPhtn)
        chkPhtn0 = find(chkPhtn,1,'first');
        idx0 = (chkPhtn0):min((chkPhtn0+maxPhtn), maxSrch);
        
        % truncate before next detection
        foo = find(recorded0(idx0),1, 'first') - 1;
        if ~isempty(foo)
            idx0 = idx0(1:foo);
        end

        % time since last event
        dt1 = cumsum(dt0(idx0));
        
        % all events fall inside of both channel deadtimes
        if all(dt1 <  t_dead_digital_min)
            chkPhtn(idx0) = 0;
            %fprintf('no event\n')
            continue
        else           
            % deremine if odd or even event
            % ** random selection of start chann should really be implimented for each pulse **
            parity=mod(sum(recorded0(1:(idx0(1)-1))),2)+1;
          
            foo = find(dt1 > t_dead_digital(i, parity),1,'first');
            if ~isempty(foo)
                %fprintf('event\n')
                recorded0(idx0(foo)) = true;
                chkPhtn(idx0(1:foo)) = false;
            else
                %fprintf('no event\n')
                chkPhtn(idx0) = 0;
            end
        end
    end
    
    foo = find(idxA);
    idxA(foo(~recorded0)) = false;
    recorded(idx) = idxA;
end

if exist('gain_params','var')
    % gain_params must be a struct with fields:
    % t
    % dt
    % n_pulses
    % 
    Nrows=length(gain_params.t);
    Ncols=gain_params.n_pulses*N_det;
    
    N_analog=(t_dead_analog/gain_params.dt);
    if N_analog>=0.25
        N_analog=ceil(N_analog);
        K_analog=[zeros(N_analog,1); ones(N_analog,1)];
    else
        K_analog=0;
    end
    N_dig=ceil(mean(t_dead_digital(:))/gain_params.dt);
    K_dig=[zeros(N_dig-1,1); ones(N_dig,1)];
    
    col=(pulse-gain_params.pulse_0)*N_det+det;
    row=ceil((t-gain_params.t(1))/gain_params.dt)+1;
    good=col>0 & col < Ncols& row > 0 & row < Nrows;
    col=col(good);
    row=row(good);
    det=det(good);
    rec_sub=recorded(good);
        
    dead_analog=conv2(full(sparse(row, col, ones(size(row)), Nrows, Ncols)), K_analog,'same')~=0;
    row=row(rec_sub); col=col(rec_sub); det=det(rec_sub);
    dead_dig= conv2(full(sparse(row, col, ones(size(row)), Nrows, Ncols)), K_dig,'same')~=0;
      
    dead_all=dead_analog | dead_dig;
 
end



