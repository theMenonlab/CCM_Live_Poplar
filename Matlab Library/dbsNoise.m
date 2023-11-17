%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve an discrete inverse problem using Direct-binary-search method
%
% Author: Ganghun Kim, Peng Wang
% Email: ganghun.kim@gmail.com
% Affiliation:  Lab for Optical Nanotechnologies,
%               Department of Electrical and Computer Engineering,
%               University of Utah,
%               Salt Lake City, UT 84112. U.S.A.
% PI: Dr.Rajesh Menon
% Copyright @ June 23,2015
%
% This is a program to reconstruct the fluorescence microscope ojbect from
% transmitted image through dielectric cylinder.
%
% Example: [x] = dbs(A, y, dI, init_type)
%
%               A           - PSF matrix
%               y           - measurement vector
%               dI          - unit perturbation of object intensity at each pixel
%               init_type   - initial solution 
%               x           - reconstructed object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x imageNos2D] = dbs(A, y, dI, init_type, imgFOV, objFOVInd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        INITIALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter sanity check
    [imageSize, objectSize] = size(A);
    if (size(y, 1) ~= imageSize)
        fprintf('Wrong dimensions! \n');
        return
    end

    % Reconstruction parameters
    dN = dI/10;
    DEBUG = false;
    method = 1;
    iterMax = 10000;                                                        % maximum number of iterations
    errorMinDelta = 1e-004;                                                % minimum error improvement
    pixelMin = 1;                                                          % minimum pixel change
    upperI = 100;
    lowerI = -0;
    imgFOVInd = find(imgFOV==1);

    % Initial solution
    %x0 = x0generator(init_type, objectSize);
    x0 = zeros(1,objectSize);

    % Reconstruct image based on initial guess
    imageOpt = zeros(imageSize, 1);
    imageNos = zeros(imageSize, 1);                                        % noise factor
    for ii = 1 : objectSize
        imageOpt = imageOpt + x0(ii) .* A(:, ii);
    end
    objectOpt = x0;                                                        % current object
    
    errorVector(1) = norm((imageOpt+imageNos)-y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                         RECONSTRUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng('default');
    for iter = 1 : iterMax

        rng('shuffle');
        sequence = randperm(objectSize + imageSize);
        %imgSequence = randperm(imageSize);
        oldError = errorVector(iter);
        objCount = 0;
        nosCount = 0;

        % Binary search
        if (method == 1)
            for ss = 1 : objectSize
                is_better = 0;
                indexRaw = sequence(ss);
                
                
                if (objectSize >= indexRaw)
                    %%
                    index = indexRaw;
                    oldI = objectOpt(index);
                    objectOpt(index) = oldI + dI;
                    
                    if (objectOpt(index) <= upperI)
                        imageNew = imageOpt + dI .* A(:, index);
                        newError = norm((imageNew + imageNos) - y);
                        if (newError < oldError)
                            imageOpt = imageNew;
                            oldError = newError;
                            is_better = 1;
                        else
                            objectOpt(index) = oldI;
                            is_better = 0;
                        end
                    else
                        objectOpt(index) = oldI;
                        is_better = 0;
                    end
                
                    if (is_better == 0)
                        objectOpt(index) = oldI - dI;
                        if (objectOpt(index) >= lowerI)
                            imageNew = imageOpt - dI .* A(:, index);
                            newError = norm((imageNew + imageNos) - y);
                            if (newError < oldError)
                                imageOpt = imageNew;
                                oldError = newError;
                            else
                                objectOpt(index) = oldI;
                            end
                        else
                            objectOpt(index) = oldI;
                        end
                    end
                
                    if (objectOpt(index) ~= oldI)
                        objCount = objCount + 1;
                    end
                    
                else
                    %%
                    if iter > -1
                        index = indexRaw - objectSize;
                        oldI = imageNos(index);
                        imageNos(index) = oldI + dN;
                        %if (imageNos(index) <= upperI)
                            newError = norm((imageOpt + imageNos) - y);
                            if (newError < oldError)
                                oldError = newError;
                                is_better = 1;
                            else
                                imageNos(index) = oldI;
                                is_better = 0;
                            end
                        %else
                        %    imageNos(index) = oldI;
                        %    is_better = 0;
                        %end

                        if (is_better == 0)
                            imageNos(index) = oldI - dN;
                            if (imageNos(index) >= lowerI)
                                newError = norm((imageOpt + imageNos) - y);
                                if (newError < oldError)
                                    oldError = newError;
                                    is_better = 1;
                                else
                                    imageNos(index) = oldI;
                                end
                            else
                                imageNos(index) = oldI;
                            end
                        end

                        if (imageNos(index) ~= oldI)
                            nosCount = nosCount + 1;
                        end
                    end
                end
                    
                    
            end
            errorMin = oldError;

        else
            disp('Wrong method! \n');
            return
        end

        % Save and output iteration results
        errorVector(iter + 1) = errorMin;
        errorChange = errorVector(iter) - errorMin;
        objPixelChange(iter + 1) = objCount;
        nosPixelChange(iter + 1) = nosCount;

            fprintf('Iteration%3d : ', iter);
            fprintf('RMS=%.4f  ', errorMin);
            fprintf('Obj_Pixel=%5d  ', objPixelChange(iter + 1));
            fprintf('Nos_Pixel=%5d  ', nosPixelChange(iter + 1));
            fprintf('\n');
            %toc;
            %figure(1); plot(objectOpt);
                    
         if (DEBUG == true)
            psfSize = 111;
            objRaw = zeros(1,psfSize^2);
            objRaw(objFOVInd) = objectOpt;
            objRaw = reshape(objRaw,[psfSize psfSize]);
            objRaw = abs(objRaw);
            figure(1); imagesc(objRaw); axis square; colormap gray; colorbar;
            
            imgXSize = size(imgFOV,1);
            imageNos2D = zeros(1,imgXSize^2);
            imageNos2D(imgFOVInd) = imageNos;
            imageNos2D = reshape(imageNos2D,[imgXSize imgXSize]);
            figure(2); imagesc(imageNos2D); axis square; colormap gray; colorbar;
        end

        % Termination conditions
        if (objPixelChange(iter + 1) < pixelMin)
            %fprintf('~~~~~ Minimum Pixel Chanage Reached ~~~~~\n');
            break;
        elseif (errorChange < errorMinDelta/1e4)
            %fprintf('~~~~~ Minimum RMS Chanage Reached ~~~~~\n');
            break;
        end
    end

    x = objectOpt';
    errorVector = errorVector(1:iter+1);
    objPixelChange = objPixelChange(1:iter+1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     RECONSTRUCTION RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (DEBUG)
        % Object location found
        figure(1); plot(x);
        title('Reconstructed x', 'FontName', 'Arial', 'FontSize', 16);

        % Capture image
        figure(2); plot(y);
        title('Measurement y', 'FontName', 'Arial', 'FontSize', 16);

        % Reconstructed image
        figure(3); plot(imageOpt);
        title('Reconstructed y', 'FontName', 'Arial', 'FontSize', 16);

        % Error of reconstruction
        figure(4); plot(imageOpt - y);
        title('Reconstruction Error', 'FontName', 'Arial', 'FontSize', 16);

        % Evolution of DBS algorithm
        figure(5);
        subplot(2,1,1);
        semilogy((0:iter), errorVector, 'ko-'); hold on;
        xlabel('Times of Iteration', 'FontName', 'Arial', 'FontSize', 14);
        ylabel('RMS Error', 'FontName', 'Arial', 'FontSize', 14);
        set(gca, 'FontName', 'Arial', 'FontSize', 12);
        subplot(2,1,2);
        semilogy((1:iter), objPixelChange(2:iter+1), 'ko-'); hold on;
        xlabel('Times of Iteration', 'FontName', 'Arial', 'FontSize', 14);
        ylabel('Pixels', 'FontName', 'Arial', 'FontSize', 14);
        set(gca, 'FontName', 'Arial', 'FontSize', 12); 
    end
    
            imgXSize = size(imgFOV,1);
            imageNos2D = zeros(1,imgXSize^2);
            imageNos2D(imgFOVInd) = imageNos;
            imageNos2D = reshape(imageNos2D,[imgXSize imgXSize]);

return
