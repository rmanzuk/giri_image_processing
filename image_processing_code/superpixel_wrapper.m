% wrapper to go through, make superpixels for all images and
%
% R. A. Manzuk Thursday, January 26, 2023 at 5:14:52 PM
% Last edited: Thursday, January 26, 2023 at 5:15:07 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% paths to all channels
in_470 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/470nm';
in_505 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/505nm';
in_530 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/530nm';
in_590 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/590nm';
in_625 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/625nm';
in_760 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/760nm';
in_940 = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/940nm';
in_365 = '/Volumes/Extreme SSD/sms_images/original_tifs/365nm';
%% make superpixels and save the label matrices
% we'll just use RGB for this
rgb_dirs = {in_625, in_530, in_470};

% number of superpixels
n_superpix = 10000;

% and make the superpixels
superpix_output_dir = '/Volumes/Extreme SSD/sms_images/superpixel_inds';
save_superpixels(rgb_dirs, superpix_output_dir, '.tif', n_superpix);
