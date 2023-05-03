## Registration

Rigid and affine registration based on tissue segmentation for creation of atlases.

## Workflow

1. Select a sample of the group to be the Objective of the Registration.
2. Using the volume of segmented tissues (preferably with manual corrections), place the 3 neck landmarks, 5 face landmarks for rough registration, and the 38 landmarks for Geometric Morphometrics.
3. Execute Script_01 to place the sample in the base of a Reference space and leave only the head of the sample. List associated volumes (cells, proliferating cells, etc) to be moved too.
4. Place same landmarks as in step 2 in the remaining samples of the group.
5. List the remaining samples in Script_02 and execute the script.
6. With all the tissue volumes rigid registered to the Objective, the affine Groupwise Registration can be performed using Script_03. The execution of the script can take time.
7. Once finished, list the all the samples of the group in Script_04 and execute the script to obtain the required mean transformation for each sample to build the tissue atlas of the group.
8. Run Script_05 for the rigid-registered associated volumes (output step 3) to apply the mean transformation.
9. Using MATLAB, in Script_06 list manually the affine registered volumes or use the command dir to find them. Run the script to generate the tissue atlas of the group.

Note: use the CSV to/from PTS to convert landmark files.
