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

function [x] = dbs(A, y, dI, init_type)
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
    DEBUG = false;
    method = 1;
    iterMax = 1000;                                                        % maximum number of iterations
    errorMinDelta = 1e-004;                                                % minimum error improvement
    pixelMin = 1;                                                          % minimum pixel change
    upperI = 100;
    lowerI = -0;

    % Initial solution
    %x0 = x0generator(init_type, objectSize);
    x0 = zeros(1,objectSize);

    % Reconstruct image based on initial guess
    imageOpt = zeros(imageSize, 1);
    for ii = 1 : objectSize
        imageOpt = imageOpt + x0(ii) .* A(:, ii);
    end
    objectOpt = x0;                                                        % current object
    errorVector(1) = norm(imageOpt-y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                         RECONSTRUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng('default');
    for iter = 1 : iterMax

        rng('shuffle');
        sequence = randperm(objectSize);
        oldError = errorVector(iter);
        count = 0;

        % Binary search
        if (method == 1)
            for ss = 1 : objectSize
                is_better = 0;
                index = sequence(ss);

                oldI = objectOpt(index);
                objectOpt(index) = oldI + dI;

                if (objectOpt(index) <= upperI)
                    imageNew = imageOpt + dI .* A(:, index);
                    newError = norm(imageNew - y);
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
                        newError = norm(imageNew - y);
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
                    count = count + 1;
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
        pixelChange(iter + 1) = count;

        %if (DEBUG == true)
            fprintf('Iteration%3d : ', iter);
            fprintf('RMS=%.4f  ', errorMin);
            fprintf('Pixel=%5d  ', pixelChange(iter + 1));    
            fprintf('\n');
            %toc;
            %figure(1); plot(objectOpt);
            %figure(1); imagesc(reshape(objectOpt,[sqrt(objectSize) sqrt(objectSize)]));
        %end

        % Termination conditions
        if (pixelChange(iter + 1) < pixelMin)
            %fprintf('~~~~~ Minimum Pixel Chanage Reached ~~~~~\n');
            break;
        elseif (errorChange < errorMinDelta/1e4)
            %fprintf('~~~~~ Minimum RMS Chanage Reached ~~~~~\n');
            break;
        end
    end

    x = objectOpt';
    errorVector = errorVector(1:iter+1);
    pixelChange = pixelChange(1:iter+1);

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
        semilogy((1:iter), pixelChange(2:iter+1), 'ko-'); hold on;
        xlabel('Times of Iteration', 'FontName', 'Arial', 'FontSize', 14);
        ylabel('Pixels', 'FontName', 'Arial', 'FontSize', 14);
        set(gca, 'FontName', 'Arial', 'FontSize', 12); 
    end

return
