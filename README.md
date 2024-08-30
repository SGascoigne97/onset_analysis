# Incomplete resection of the icEEG seizure onset zone is not associated with post-surgical outcomes
Gascoigne et al. (brief communication, accepted July 2024).
This code performs downstream analysis, comparing proportion of onset resected, onset size, and onset as a proportion of implantation across outcome groups. 
Within this folder **run_all_analyses.m** calls all other scripts and will produce tables and figures which were used in this paper. This script runs analyses for all parcellations (channel-wise, Lausanne-120, Lausanne-250) - changing line 45 to:
  for parc = "roi_120"
In this case the comparison across parcellations is not possible (i.e. code from line 311 cannot be run)

In the main analysis, we looked at onsets based on the Lausanne-120 parcellation scheme, this code can be easily adapted to produce results using channel-wise onsets or the Lausane-250 parcellation scheme. 

#### Setting parameters for analysis
- *parc*: The parcellation scheme of interest (can be set to "chan", "roi_120", or "roi_250"). You only need to set this parameter if you are not iterating through parcellation schemes (i.e., lines 40 and 86 remain commented out). The code is currently set to perform analysis using the Lausanne-120 parcellation scheme (Line 17: *parc = "roi_120";*)
- *n_perm*: The number of permutations to perform to create a "chance" distribution of AUCs from which we can compute the p-value. The smaller this number is set, the faster the code will run, however we suggest that this should not be set below 500 if the p-values are to be reliable. Our work used 1000 permutations to build the distribution of AUCs (Line 18: *n_perm = 1000;*)
- *out_thresh*: The maximum ILAE values that is labelled as a 'favourable' surgical outcome. In this work we have set this value to 2 such that ILAE 1-2 and ILAE 3+ are labelled as having 'favourable' and 'unfavourable' surgical outcomes, respectively (Line 20: *out_thresh = 2;*)
- *consensus_thresh*: The threshold for inclusion of regions in the consensus onset (i.e. one onset per individual which captures regions that are prevalent across seizures). Here we will include regions present in at least half of the subject's seizures (Line 22: *consensus_thresh = 0.5;*)
- *most_resec_thresh*: The threshold above which subjects are considered as having 'most' of their onset resected. Here we have set this as 50%, so any resections removing at least half of the onset is labelled as 'most resected' (Line 25: *most_resec_thresh = 0.5;*)
- *det_meths*: This selects which methods of onset detection will be considered (can be set to "clo", "imprint", or "EI"). In our main text, we used clinically labelled onsets (CLO) and automatically labelled onsets (ALO) based on the imprint algorithm. Supplementary analyses additionally included automatically detected onsets using the Epileptigenicity Index (Bartolomei et al., 2008) (Line 28:*det_meths* = ["clo", "imprint", "EI"];*)

#### Saving figures
The code saves figures for each step of the analysis, the user can determine if they want the figures to be saved (save_fig = 1) or not (save_fig = 0) on line 37. 
The location where the figures are stored can be changed, as this code produces the visualisations used in *Figure Two*, we have set the save location accordingly (Line 34: *save_folder = '../figures/paper_figures/Figure 2/'*).
Additionally, it is possible to change the filetype used for each of the figures (Line 37: *file_type = "svg";*), see MATLAB documentation to see the filetypes available). 

#### Step-by-step downstream analysis
##### 3.1 icEEG seizure onset regions tend to be resected, but more complete resections are not associated with more favourable surgical outcomes
**is_onset_resected** The proportion of the seizure onset zone (CLO and ALO) that was subsequently resected was computed and compared across outcome groups.  We elected to omit the volume calculations from this work as the volume of regions may not accurately represent the true volume of the onsets and resections, the code has been left in this repository for completeness. Output figures will include comparisons using both the number of regions and the estimated volume of onsets.

##### 3.2 Larger onsets are not associated with surgical outcomes
**onset_resec_as_prop**
Here we compute the size of seizure onsets (both CLO and ALO) as a count of regions and proportion of implantation and compare across surgical outcome groups. Again, we elected to omit the volume calculations from this work as the volume of regions may not accurately represent the true volume of the onsets and resections, the code has been left in this repository for completeness. Output figures will include comparisons using both the number of regions and the estimated volume of onsets.

#### Figures
*Figure Two* is based on the figures generated when running **run_all_analyses.m**. The panels used in Figure 2 are saved in the following locations:
- A: *figures/paper_figures/Figure 2/onset_resected/roi_120_Perc_violin.svg*
- B: *figures/paper_figures/Figure 2/onset_size/roi_120_violin_count.svg*
- C: *figures/paper_figures/Figure 2/resection_size/roi_120_violins.svg*


