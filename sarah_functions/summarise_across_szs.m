function [arr] = summarise_across_szs(mat, min_sz)
% Take matrix of automatically detected onsets, and summarise into array
% using a threshold for the minimum number of seizures in which a region is
% detcted as onset

% Sarah J Gascoigne 05/06/2023

% input:
%   - mat: binary matrix of onsets across seizures
%   - min_sz: numeric denoting the minimum number of seizures required for
%   region to be included
%   

% output
%   - arr: binary array of regions in onset surpassing threshold

    arguments
        mat (:,:) logical % binary matrix of onsets across seizures
        min_sz (1,1) {mustBeNumeric(min_sz)} % numeric denoting the minimum number of 
                             % seizures required for region to be included
       
    end

    if(size(mat,2)==1)
        fprintf("Only one seizure recorded, summary is not required \n")
         arr = mat;
        return 
    end

    arr = sum(mat,2)>=min_sz;
end