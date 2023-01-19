function [adjustments] = mean_to_lcc(mean_im)
% This function takes a mean image from a dataset and fits a surface to it
% to best approximate lcc adjustments. The function asks the user to click
% the outline of a polygon that is the best-constrained part of the mean
% image for masking in order to avoid using misrepresentative areas like
% those that can occur in the corners of some images.
%
% IN:
% mean_im: 2d double matrix with the mean image you would like to fit a
% surface to for lcc correction
%
% OUT:
% adjustments: 2d double matrix of equal size to the input image with the
% needed adjustments to achieve even lighting
%
% Written by R.A. Manzuk
% 10/18/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % we'll do this on a smoothed and downsampled version of the mean image,
    % for efficiency
    downsamped_mean = imresize(mean_im, 0.05);
    smoothed_mean = imgaussfilt(downsamped_mean,size(downsamped_mean,1)/10);
    
    % define a mesh grid for the image
    x = 0:size(smoothed_mean,2)-1;
    y = 0:size(smoothed_mean,1)-1;
    [X,Y] = meshgrid(x,y); 
    
    % have the user input a polygon where the mean image is roughly well
    % defined
    % set up the empty outer polygon coordinates
    outer_poly = [];
    
    disp('please click the outline of a rough polygon that contains the well-defined area.')
    imagesc(smoothed_mean)
    colorbar
    ax = gca;
    ax.Toolbar.Visible = 'off';
    hold on
    
    title('Polygon tracing'); 
    y_coords = [];
    x_coords = [];
    n = 0;
    while true
        [x_i,y_i] = ginput(1);
        if isempty(x_i) ; break; end
        n = n+1;
        x_coords(n) = x_i(1);
        y_coords(n) = y_i(1);
        plot(x_coords,y_coords,'r','LineWidth',1)
        drawnow
    end
    
    % make this polygon into a mask for the mean image
    mask = poly2mask(round(x_coords),round(y_coords), size(smoothed_mean,1), size(smoothed_mean,2));
    
    % to do the surface fitting, etc, we'll need vector versions of everything
    x_vec = X(:);
    y_vec = Y(:);
    mean_vec = smoothed_mean(:);
    mask_vec = mask(:);
    
    % fit the surface, we'll assume 3rd order polynomials in both axes
    sf = fit([x_vec(mask_vec),y_vec(mask_vec)],mean_vec(mask_vec),'poly22');
    
    % and evaluate it over all x,y....not just the max
    fit_mean = sf(x_vec,y_vec);
    
    % let's scale our solved surface such that the max value is 1
    scaled_mean = fit_mean + (1-max(fit_mean));
    
    % now just have to reshape and resize
    mean_reshaped = reshape(scaled_mean,size(smoothed_mean,1),size(smoothed_mean,2));
    adjustments = imresize(mean_reshaped,[size(mean_im,1),size(mean_im,2)]);
end
