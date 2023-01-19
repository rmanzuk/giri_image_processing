function [var_pct_explained, max_var_inds] = assess_anisotropy(response_cols, orient_inds, n_filters)
% Knowing that the percentage of variance explained by the
% maximally-aligned filter is a good metric for image anisotropy (based
% upon testing with icnofabric index images), we can apply the metric to
% an image.
%
% INPUT 
%
% response_cols: matrix with the normalized responses of a filter bank on
% an image. Each colum represents a different filter, and each row is a
% different pixel.
%
% orient_inds: not all filters in the filter bank may be oriented (some
% could be blobs, etc.). This is a logical vector with the length of the
% number of filters in the filter bank. True indicates an oriented filter.
%
% n_filters: number of filters to be considered in the max variance
% calculation. 6 usually is about right
%
% OUTPUT
%
% var_pct_explained: Percentage of varience explained by the input
% maximally varient filters.
%
% max_var_inds: indices of the maximally variant filters in the filter
% bank. In order of variance.
% 
% Written by R.A. Manzuk
% 11/29/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % grab the oriented filter responses
    oriented_responses = response_cols(:,orient_inds);
        
    % and get the filters that have the highest variance
    [max_vars,max_var_inds] = maxk(var(oriented_responses,'omitnan'),n_filters);
    
    % what percentage of the variance in filter responses is accounted for
    % by just those?
    var_pct_explained = 100 * sum(max_vars)/sum(var(response_cols,'omitnan'));
    
end