function [img_hdr] = getImageHDR(path, Exp, filename)
    
    load HDR_Orca;
    
    jmax = size(Exp,2);
    for j = 1:jmax
        if strcmp(filename, 'rawImage')
            temp = path{j};
        elseif strcmp(filename, 'exp')
            temp = imread([path,'Exp',num2str(Exp(j)),'.tif']);
        elseif strcmp(filename, 'single')
            temp = imread([path]);
        elseif strcmp(filename, 'legacyCal')
            temp = imread([path,'Exp',num2str(Exp(j)),'/',filename]);
        else
            temp = imread([path,filename,'_Exp',num2str(Exp(j)),'.tif']);
        end
        if j == 1
            img_hdr = zeros(size(temp));
            normalizer = zeros(size(temp));
        end
        temp = temp - darkLevel;
        temp(temp<0) = 0;
        temp(temp>whiteLevel) = whiteLevel;
        Z = temp+1;
        B(j) = log(Exp(j));
  
        if size(img_hdr) ~= size(w(Z))
            error('HDR Error. Size : [%d %d], Cell size [%d %d]: .',size(img_hdr), size(path));
        end
        
        img_hdr = img_hdr + w(Z).*(g(Z) - B(j));
        normalizer = normalizer + w(Z);
    end
    
    img_hdr  = img_hdr ./ normalizer;
    img_hdr  = exp(img_hdr);
    img_hdr(isnan(img_hdr)) = 0;
    
    clear relExp img;
    clear w z B Z g i j k jmax normalizer path temp filename DC;
    
end