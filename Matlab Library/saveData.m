%% Variable Needed
% Every variables to be saved
% savefile : path to save the file

if ~saveObjPro
    for f = 1:size(frames,2)
        frames(f).objPro = 0;
    end
end

if ~exist('calID')
    calID = 0;
    calPath = 'NA';
end

save(filename,'frames','calPath','calID','saveObjPro','-v7.3');
