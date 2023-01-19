function [parsed_exif] = get_parse_exif(path_to_im) 
% This function will use the exiftool to extract exif data for an image
% file, parse the raw character version of the data, and returns a
% structure with the organized exif data
% 
% This function requires exif tool: https://oliverbetz.de/pages/Artikel/ExifTool-for-Windows
%
% IN:
% path_to_im: string with the full path to the image
%
% OUT:
% parsed_exif: structure with all exif fields and data
%
% Written by R.A. Manzuk 
% 10/06/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % just start by getting the exif
    exif_raw = getexif(path_to_im);
    % make an empty structure to hold the parsed exif
    parsed_exif = struct; 
    % split the character array into lines
    lines = split(exif_raw,char(10));

    % iterate through to place all fields/data in the structure
    for j = 1:length(lines)-1
       % split the lines by the colon
       parts = split(lines{j}, ':');
       % If there are just two parts, enter straight away
       if length(parts) == 2
           parsed_exif.(deblank(parts{1})) = deblank(parts{2});
       % if more than two parts, join the later parts
       elseif length(parts) > 2
           parsed_exif.(deblank(parts{1})) = deblank(strjoin(parts(2:length(parts)), ':')); 
       end
    end
end