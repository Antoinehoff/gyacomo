function [field, Ts2D] = compile_results_2Da(DIRECTORY,JOBNUMMIN,JOBNUMMAX,fieldname)
CONTINUE = 1;
JOBNUM   = JOBNUMMIN; JOBFOUND = 0;
% field
field    = [];
Ts2D    = [];
ii = 1;
while(CONTINUE)
    filename = sprintf([DIRECTORY,'outputs_%.2d.h5'],JOBNUM);
    % Check presence and jobnummax
    if (exist(filename, 'file') == 2 && JOBNUM <= JOBNUMMAX)
        %test if it is corrupted or currently running
        try
            openable = ~isempty(h5read(filename,'/data/var2d/time'));
        catch
            openable = 0;
        end
        if openable
        % load field %%
        [ F, T, ~] = load_2Da_data(filename, fieldname);
        field  = cat(4,field,F);
        Ts2D  = cat(1,Ts2D,T);
        ii = ii + 1;
        JOBFOUND = JOBFOUND + 1;
        end
    elseif (JOBNUM > JOBNUMMAX)
        CONTINUE = 0;
    end
    JOBNUM   = JOBNUM + 1;
end

if(JOBFOUND == 0)
    disp('no results found, please verify the paths');
    return;
end

end