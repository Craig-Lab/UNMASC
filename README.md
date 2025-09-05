Code to run the model described in Bistodeau-Gagnon, Craig, and Michaels. 

Each folder contains the code necessary to reproduce all the figures. Each folder can be run independently of the other. 

**DATA**  
The data analyzed in this study were obtained from a publicly available source with DOI Michaels (10.1126/sciadv.abn5522, Tables S1-S7). The Data folder contains all the cytokine concentrations in ng/mL and the polynomial coefficients for the dose-responses from Michaels et al. (2021).

**INSTRUCTIONS**  
In the code, it is worth mentioning that we refer to D7 as the T cell differentiation stage and D14 as the T cell maturation stage, as per Michaels et al., 2021 (10.1126/sciadv.abn5522). 
Cells types are denoted as follows:
	PROT: ProT cells  
	CD4ISP: CD4 single positive cells  
	CD3N: Double positive CD3- cells  
	CD3P: Double positive CD3+ cells  
	SP: CD8 single positive cells  

In each folder for figures and supplementary figures, run file containing "commands" in the title to execute the code. 
When necessary, a folder “Results” or “Figures” is included in the folder to store results files (.m files) and figures (in PDF). 

**SOFTWARE**  
All of the code was written using MATLAB R2024b. The optimization toolbox is required to perform the fitting. 
