# Cas-integrin-paper-2023-Cooper-lab
 Quantify dynamics of protein accumulation at adhesion clusters

__Cas-integrin-paper-2023-Cooper-lab__ contains MATLAB scripts used to quantify the time shift between red (CasmSc) and green (YFP-VCL, b1-EctopH, ITGB1GFP or CrkmGL) cluster formation  from TIRF microscopy datasets described in Saurav *et al*. The pipeline involves the following steps: drift correction, image preprocessing and denoising, cluster segmentation and tracking, intensities extraction and normalization, and signal analysis.

__*adhesiondynamicstirf.m*__ measures the time shift in the accumulation of two different proteins at the same adhesion cluster. It relies on __*bfmatlab*__, the bio-format package for MATLAB that can be freely downloaded at:

https://www.openmicroscopy.org/bio-formats/downloads/

*__ImregStack.m__* corrects the image stack for stage-drifting using translation transformation.

*__progressText.m__* is a helper function developped by the Danuser lab to show progress of a loop as text on the screen. It is employed by both *adhesiondynamicstirf()* and *ImRegStack()*.

*__DemoFATIRF.tif__* is an example of TIRF dataset that can be used for demo purposes.
