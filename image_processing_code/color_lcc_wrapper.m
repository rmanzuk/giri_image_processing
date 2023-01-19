% wrapper to go through and lcc + color correct images
%
% R. A. Manzuk 10/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set everything up
% paths to all lcc images
lcc_470 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/470.tif';
lcc_505 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/505.tif';
lcc_530 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/530.tif';
lcc_590 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/590.tif';
lcc_625 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/625.tif';
lcc_760 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/760.tif';
lcc_940 = '/Volumes/Extreme SSD/sms_images/standards/lcc_card/940.tif';
% paths to  inputs and outputs for lcc correction
in_470 = '/Volumes/Extreme SSD/sms_images/original_tifs/470nm';
out_470 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/470nm';
in_505 = '/Volumes/Extreme SSD/sms_images/original_tifs/505nm';
out_505 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/505nm';
in_530 = '/Volumes/Extreme SSD/sms_images/original_tifs/530nm';
out_530 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/530nm';
in_590 = '/Volumes/Extreme SSD/sms_images/original_tifs/590nm';
out_590 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/590nm';
in_625 = '/Volumes/Extreme SSD/sms_images/original_tifs/625nm';
out_625 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/625nm';
in_760 = '/Volumes/Extreme SSD/sms_images/original_tifs/760nm';
out_760 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/760nm';
in_940 = '/Volumes/Extreme SSD/sms_images/original_tifs/940nm';
out_940 = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs/940nm';

%% get the mean images for lcc

mean_im_470 = calc_mean_im(in_470,'.tif');
mean_im_505 = calc_mean_im(in_505,'.tif');
mean_im_530 = calc_mean_im(in_530,'.tif');
mean_im_590 = calc_mean_im(in_590,'.tif');
mean_im_625 = calc_mean_im(in_625,'.tif');
mean_im_760 = calc_mean_im(in_760,'.tif');
mean_im_940 = calc_mean_im(in_940,'.tif');
%% use mean ims to assess all lcc and make the corrections
% assess adjustments
adjustments_470 = mean_to_lcc(mean_im_470);
adjustments_505 = mean_to_lcc(mean_im_505);
adjustments_530 = mean_to_lcc(mean_im_530);
adjustments_590 = mean_to_lcc(mean_im_590);
adjustments_625 = mean_to_lcc(mean_im_625);
adjustments_760 = mean_to_lcc(mean_im_760);
adjustments_940 = mean_to_lcc(mean_im_940);

% make adjustments 
apply_lcc(in_470, out_470, '.tif', adjustments_470);
apply_lcc(in_505, out_505, '.tif', adjustments_505);
apply_lcc(in_530, out_530, '.tif', adjustments_530);
apply_lcc(in_590, out_590, '.tif', adjustments_590);
apply_lcc(in_625, out_625, '.tif', adjustments_625);
apply_lcc(in_760, out_760, '.tif', adjustments_760);
apply_lcc(in_940, out_940, '.tif', adjustments_940);

%% assess the gray values

cc_dir = '/Volumes/Extreme SSD/sms_images/standards/color_card';
[gray_vals,wavelengths] = assess_colorbalance(cc_dir, '.tif');

bit_depth = 16;
input_dir = '/Volumes/Extreme SSD/sms_images/lcc_corrected_tifs';
output_dir = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs';
apply_colorbalance(gray_vals,wavelengths,bit_depth,input_dir,output_dir,'.tif');

%% make some rgb jpegs.
imlog = readtable('/Volumes/Extreme SSD/sms_images/image_logs/masterlog.xlsx');

red_folder = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/625nm';
green_folder = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/530nm';
blue_folder = '/Volumes/Extreme SSD/sms_images/color_corrected_tifs/470nm';
im_in_ext = '.tif';
im_out_ext = '.jpg';
name_map = [imlog.file_name, imlog.actual_sample];
destination_folder = '/Users/ryan/Dropbox (Princeton)/reef_survey_project/nevada_qgis/rgb_jpgs';
scaling = 0.15;

make_rgbs(red_folder, green_folder, blue_folder, destination_folder, im_in_ext, im_out_ext, scaling, name_map);

%% make some masks based upon areas with bad band ratios
% we'll use all the visible wavelengths for this, but not yellow be cause
% it can get weird
mask_dirs = {in_470, in_505, in_530, in_625};

% and make the masks
mask_output_dir = '/Volumes/Extreme SSD/sms_images/masks';
background_masks(mask_dirs,mask_output_dir,'.tif');