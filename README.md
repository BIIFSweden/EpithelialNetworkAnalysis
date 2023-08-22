# Structural Analysis of the Cervical Epithelial Tissue

### Overview

### Software Requirements

The software listed below should be installed before running the scripts available in this repository.

* [Fiji](https://fiji.sc)
* [MorphoLibJ plugin](https://imagej.net/plugins/morpholibj) for Fiji

* [Matlab](https://se.mathworks.com/products/matlab.html) is required to run the script that generates the *neuriteness* image. The library *vesselness2d* is used to generate the enhanced neuriteness networks and is available in this [link](https://github.com/BoguslawObara/vesselness2d).

### Usage

1. Download the Git repository for this project.
2. Navigate to the downloaded git repository directory
3. To run the Fiji scripts (*.ijm*), open Fiji and go to Plugins – Macros – Edit... and browse the corresponding file. 
4. Set the parameters accordingly
5. Execute the macro by pressing the *Run* button

### Parameter set

### *Neuriteness* script

The Matlab script *batch_neuriteness.m* generates the *neuriteness* images for each fluorescent channel of each SNI using the method proposed by: Obara, Boguslaw, et al. "Contrast-independent curvilinear structure detection in biomedical images." *IEEE Transactions on Image Processing* 21.5 (2012): 2572-2581.

To run *batch_neuriteness.m*, load the script in Matlab and then update the corresponding local path in the *addpath* command, which should be linked to the [*vesselness2d* library](https://github.com/BoguslawObara/vesselness2d). The corresponding neuriteness images of each fluorescent channel are saved in the subfolders of the input directory. 

### Support

If you find a bug, please [raise an issue](https://github.com/BIIFSweden/EpithelialNetworkAnalysis/issues/new)

### Authors

[SciLifeLab BioImage Informatics Facility (BIIF)](https://biifsweden.github.io/) project lead: Gisele Miranda

### Licence



