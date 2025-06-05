# TRansposon Assisted Chromatin Kinetic Imaging Technology (TRACK-IT)
[![DOI](https://zenodo.org/badge/939114102.svg)](https://doi.org/10.5281/zenodo.15604035)

Joo Lee<sup>1</sup>, Liang-Fu Chen<sup>1</sup>, Simon Gaudin<sup>1,2</sup>, Kavvya Gupta<sup>1</sup>, Andrew Spakowitz<sup>3</sup>, Alistair Nicol Boettiger<sup>1</sup>

<sup>1</sup>Department of Developmental Biology, Stanford University
<sup>2</sup>Department of Genetics, Stanford University
<sup>3</sup>Department of Chemical Engineering, Stanford University

Welcome.
This repository contains code related to our project on  TRansposon Assisted Chromatin Kinetic Imaging Technology (TRACK-IT).
TRACK-IT is an approach for rapidly inserting and mapping pairs of fluorescent tags in the genomes of cells for live-cell imaging, across genomic length scales and diverse chromatin states.
The purpose of this repository is to host the software we use in the processing and analyses of these data, including the links to other repositories that we build on to accomplish these analyses.
Details of the analysis protocol will be available in our manuscript (currently in preparation). 

Questions should be addressed to Alistair Boettiger, boettiger@stanford.edu

## Dependencies:

1. Spotfinding using storm-analysis - efficient spot-fitting algorithms for dense data.
https://github.com/ZhuangLab/storm-analysis/releases

2. Particle tracing using ORCA-public https://github.com/BoettigerLab/ORCA-public - contains functions called by the TRACK-IT analysis.  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15603850.svg)](https://doi.org/10.5281/zenodo.15603850)


3. Polymer simulations from polychrom https://github.com/BoettigerLab/polychrom - software for molecular dynamics simulations of loop extrusion. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7698987.svg)](https://doi.org/10.5281/zenodo.7698987)
forked from https://github.com/open2c/polychrom. 

## Organization
### Particle Tracking: 
The script Demo_ProcessTracesBatch.m provides an overview of the paticle tracking for trace assembly process. 
This script coordinates the pipeline of analysis from spot-fitting to individual cell traces. Most of the action happens inside separate functions, which are included in the ORCA-public repository.
This includes correcting for measured offsets in x,y,and z, correcting chromatic aberration, linking spots into traces while filtering background noise and accomidating missed detection events, segregating sister chromatids in G2 traces, converting data to 3D distances in nanometers, and exporting the data as tables.

### ChromaticCorrectionTool
Contains a graphical user interface App (written in the matlab App builder) to compute 3D chormatic corrections from 3D stacks of fluorescent beads embedded in 3D gels.

### FigureCode
Contains scripts used to produce the graphs and analyses presented in the figures of our manuscript

### ParameterFiles
Contains example parameter files used in mufit-analysis from the storm-analysis toolbox, as described in the methods. 

## Data availability
 Live cell tracking data is available on Zenodo: doi:10.5281/zenodo.15603650 and has been submitted for hosting at the 4DN data portal and is under review.  Accession numbers will be updated when available.
