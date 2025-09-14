Code to run the model described in Bistodeau-Gagnon, Craig, and Michaels. 

### System requirements
All of the code was written and tested using MATLAB R2024b. The following toolboxes are required:  
	- Optimization Toolbox  
 	- Statistics and Machine Learning Toolbox  
  	- Symbolic Math Toolbox  

### Installation guide
Download all the files in this repository. The typical install time is < 1 min. 

### Demo and instructions for use
Each folder contains the code necessary to reproduce all the figures. Each folder can be run independently of the other.  
In the code, it is worth mentioning that we refer to D7 as the T cell differentiation stage and D14 as the T cell maturation stage, as per Michaels et al., 2021 (10.1126/sciadv.abn5522). 
Cells types are denoted as follows:  
	PROT: ProT cells  
	CD4ISP: CD4 single positive cells  
	CD3N: Double positive CD3- cells  
	CD3P: Double positive CD3+ cells  
	SP: CD8 single positive cells  

In each folder for figures and supplementary figures, run file containing "commands" in the title to execute the code. 
When necessary, a folder “Results” or “Figures” is included in the folder to store results files (.m files) and figures (in PDF). 
The expected run time for each folder is <10 min.

### Data
The data analyzed in this study were obtained from a publicly available source (10.1126/sciadv.abn5522, Tables S1-S7). The Data folder contains all the cytokine concentrations in ng/mL and the polynomial coefficients for the dose-responses from Michaels et al. (2022).



