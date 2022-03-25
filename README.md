# Critical Role of Insertion Preference for Invasion Trajectory of Transposons

## Munasinghe, Springer, and Brandvain (in prep)

This repository contains all of the SLiM scripts used to explore how the transpostion probability and insertion preference of a TE family influence the evolution of mean TE copy number and host populaton size, allowing for extinction. A set of 5 distinct models were developed to assess the outcome of introducing a TE, which can replicate freely, into a naïve population. These models varied ether the genomic architecture of the host population or specific aspects of TE biology. 

| Model      | Number of Chromosomes | Recombination Rate | Excision | Non-Autonomous Elements |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| 1      | One       | 1.0e-08      | No       | No      |
| 2      | One       | 1.0e-08      | Yes       | No      |
| 3      | One       | 1.0e-08      | No       | Yes      |
| 4      | One       | 1.0e-08      | Yes       | Yes      |
| 5      | Five       | 1.0e-05      | No       | No      |

The folder \textbf{SLiM_Models} contains the associated SLiM scripts named following this conventon (ModelX.slim). The header of each of these files provides details about the genomic architecture and TE biology considered. 
