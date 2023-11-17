startTime = datestr(now);

startTime = datestr(now);
save([path,'/','calMetafile.mat'],'orcaExp', 'Xinit', 'Yinit', 'Zinit', 'xCount','yCount','zCount','xStep','yStep','zStep', ...
              'binning','AmCrop', 'OrcaCrop', 'AmGain','LEDPower','path','startTime', 'AmExp', 'CalRepeat', 'CalWait');
          
% clear('calExp','calAvg','xCount','yCount','zCount','xStep','yStep','zStep', ...
%               'binning','hPos','vPos','imageSize','extPower','startTime');