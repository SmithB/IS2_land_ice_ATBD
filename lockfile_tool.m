function varargout=lockfile_tool(varargin)

if nargout>0
    [varargout{1:nargout}]=feval( varargin{:});
else
    feval(varargin{:});
end

function lock_dir(file_name)
[thedir, thefile]=fileparts(file_name);

this_PID=feature('getPID');

lockdir_file=[thedir,'/LOCK_DIR_', num2str(this_PID)];
if exist(lockdir_file,'file')
    delete(lockdir_file)
end

keep_waiting=true;
while keep_waiting
    lock_list=dir([thedir,'/LOCK_DIR_*']);
    while ~isempty(lock_list)
        % there's already a lock on this directory, wait until there isn't
        pause(0.1)
        lock_list=dir([thedir,'/LOCK_DIR_*']);
    end
    % The directory looks unlocked.  write to a lock file.
    fid=fopen(lockdir_file,'w'); fwrite(fid, 1,'char'); fclose(fid);
    % check if ours is the lock file with the smallest date
    lock_list=dir([thedir,'/LOCK_DIR_*']);
    if length(lock_list) > 1
        LockDate=[lock_list.datenum]; 
        Our_file=dir(lockdir_file);
        
        if min(LockDate) < Our_file.datenum
            % another process has won the race.  Wait until it
            % is done.
            delete(lockdir_file);
        else
            % we've won the race.
            keep_waiting = false;
        end
    else
        % we're the winner
        keep_waiting=false;
    end
end

disp('');

function unlock_dir(file_name)
[thedir, thefile]=fileparts(file_name);
this_PID=feature('getPID');

lockdir_file=[thedir,'/LOCK_DIR_', num2str(this_PID)];
delete(lockdir_file);


function lock_name=lockfile_name(file_name)
lock_name=[file_name,'_LOCK'];

function status=check(file_name)
lock_name=lockfile_name(file_name);
status=exist(lock_name,'file');

function status=lock(file_name)
lock_dir(file_name);

lock_name=lockfile_name(file_name);
status=exist(lock_name,'file');

if status==0
    fid=fopen(lock_name,'w');
    fwrite(fid, 1,'double');
    fclose(fid);
end
unlock_dir(file_name);

function status=unlock(file_name);
lock_name=lockfile_name(file_name);
status=exist(lock_name,'file');

if status > 0
    delete(lock_name);
end


