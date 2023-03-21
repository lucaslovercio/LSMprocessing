# Light Sheet Microscopy image processing for developmental biology

Lucas D. Lo Vercio et al. (lucasdaniel.lovercio@ucalgary.ca)

Cumming School of Medicine, University of Calgary (Calgary, AB, Canada)

## Introduction

It is provided:
* Dataset of annotated DAPI-stained and PHH3-stained slices.
* Scripts and pretrained U-net for segmentating mesenchyme and neural ectodem in Lightsheet Microscopy (LSM)  DAPI-stained images.
* Pretrained network for segmenting cells in DAPI-stained slices using the Fiji-plugin by Falk et al. (2019).
* Pretrained network for segmenting proliferating cells in PHH3-stained slices using the Fiji-plugin by Falk et al. (2019).

## Datasets

The used datasets and trained U-nets (H5 files) can be found in:

https://cumming.ucalgary.ca/lab/morpho/code-and-data

## Software requirements

* Tissue segmentation: Python 3, TensorFlow, Keras, imgaug (https://github.com/aleju/imgaug).
* Cells and proliferating cells segmentation: ImageJ-Fiji, and Fiji-plugin by Falk et al. (2019) installed and running.
* For tiling and merging: MATLAB. This can be replaced by your own scripts.

## Workflow

1. Tile the DAPI-stained slices with defined overlapping. You can use the provided MATLAB_Scripts/Script_01_scriptTileVolumeNB.m
2. Segment the tissues using our proposed U-net (Python_Scripts/Script_03_scriptSegmentTiles.py)
3. Merge the resulting tiles with tissues segmented. You can use the provided MATLAB_Scripts/Script_04_scriptMergeProcessedTilesNB.m
4. Segment the cells in the DAPI-stained tiles generated in step 1, using the corresponding trained architecture and the script provided in ImageJ_Fiji_Scripts/segment_folder.ijm
5. Merge the resulting tiles with cells segmented. You can use the provided MATLAB_Scripts/Script_04_scriptMergeProcessedTilesNB.m
6. Tile the PHH3-stained slices with defined overlapping. You can use the provided MATLAB_Scripts/Script_01_scriptTileVolumeNB.m
7. Segment the cells in the PHH3-stained tiles, using the corresponding trained architecture and the script provided in ImageJ_Fiji_Scripts/segment_folder.ijm
8. Merge the resulting tiles with proliferating cells segmented. You can use the provided MATLAB_Scripts/Script_04_scriptMergeProcessedTilesNB.m
9. (Optional) Create volumes as multipage TIFF using MATLAB_Scripts/Script_05_scriptCreateVolumeFromSlices.m
10. (Optional) If you want to digitally remove parts of the volume, we suggest to use the tissue segmentation volume for this task and then mask the cell and proliferative cells segmentation volumes. For cropping, you can use 3D Slicer (https://www.slicer.org/).

## License

MIT License

Copyright 2020 Lucas D. Lo Vercio

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### Citing

L. D. Lo Vercio et al., "Segmentation of Tissues and Proliferating Cells in Light-Sheet Microscopy Images of Mouse Embryos Using Convolutional Neural Networks," in IEEE Access, vol. 10, pp. 105084-105100, 2022, doi: 10.1109/ACCESS.2022.3210542.

...

@ARTICLE{LoVercio2022,
  author={{Lo Vercio}, Lucas D. and Green, Rebecca M. and Robertson, Samuel and Guo, Sienna and Dauter, Andreas and Marchini, Marta and Vidal-García, Marta and Zhao, Xiang and Mahika, Anandita and Marcucio, Ralph S. and Hallgrímsson, Benedikt and Forkert, Nils D.},
  journal={IEEE Access}, 
  title={Segmentation of Tissues and Proliferating Cells in Light-Sheet Microscopy Images of Mouse Embryos Using Convolutional Neural Networks}, 
  year={2022},
  volume={10},
  number={},
  pages={105084-105100}
}

...


### References

Falk, T. et al., U-net: deep learn-ing for cell counting, detection, and morphometry, Nature methods (2019).
