# An ICR Tutorial
This respository contains all the files and scripts necessary to participate in the FTICR-MS data analysis demonstration at the SBR SFA Community Workshop 2019. 
Moreover, this repository will hopefully provide insight into FTICR-MS analysis more broadly.

Here is a brief outline of the included files:

<b>Carbon_Efficiency.R:</b> This script will calculate carbon-use efficiency for each carbon content in a dataset and will generate one figure.<br>
<b>FTICR Analysis Walkthrough.docx:</b> This Word document has step-by-step instructions which will walk through how to analyze the provided FTICR data<br>
<b>FTMS_Analysis.R:</b> This script is a necessary pre-processing step for a provided FTICR report in order for other scripts to function. This will remove isotopic peaks and limit the mass range of peaks to 200-900 m/z. This script is based upon the <a href=https://github.com/EMSL-Computing/ftmsRanalysis/tree/master>'ftmsRanalysis' R package</a><br>
<b>HJ_Andrews_Tutorial_Report.csv:</b> An example FTICR report generate from Watershed 1 within the HJ Andrews Experimental Forest, Oregon<br>
<b>SBR_SFA_ICR_Tutorial.pptx:</b> The PowerPoint shown at the Watershed Workshop ICR Tutorial<br>
<b>Transformation_Analysis.R:</b> This script will perform a "transformation analysis" which will identify putative biochemical reactions within a dataset given a transformation database. 
Such analyses are outlined in <a href=https://doi.org/10.1007/s11306-006-0029-z>Breitling et al., 2006 - Metabolomics</a> and <a href=https://doi.org/10.1002/2017JG003967>Graham et al., 2017 - JGR: Biogeosciences</a>. This script will generate two figures.<br>
<b>Transformation_List_Amino_Acid.csv:</b> The transformation database for the transformation analysis; only identifies putative amino acid transformations.<br>
<b>install_packages_commands.txt:</b> A file to allow copying-and-pasting of R install package commands.<br>
