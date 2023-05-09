# Fontaine Rynearson 2023 (in review at L&O)
Associated code and custom database used for Fontaine and Rynearson 2023 manuscript 

Files:
1. read_threshold.R: Takes the mock community counts and taxonomy files and determines an appropriate read cut-off threshold in terms of read % abundance. The threshold we determined for this manuscript is 0.075%

2. Script_for_Fontaine_Rynearson.R: Main script used to process all data and create figures for this manuscript. The figures are commented out based on figure name in manuscript

3. Skel_Thal_custom_DB.fa: The custom database used to assigned taxonomy to Skeletonema and Thalssiosira ASVs that were unassigned at the species level. This custom database includes the full 18S sequences when possible and was built from previous literature (Canesi and Rynearson, 2016 and Rynearson et al. 2020). Taxonomic assignment only included the use of the V4 region of these sequences.


