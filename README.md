# An FTICR-MS Data Analysis Tutorial
This respository contains all the files and scripts necessary to participate in various FTICR-MS data analysis demonstrations. Moreover, this repository should provide some insight into FTICR-MS analysis more broadly. This will not outline all possible routes of analysis, but only demonstrate only a few. Many of these scripts were tailored for the tutorial specified, so they might not be completely compatible with other datasets without some additional modification. If there are any questions, please refer them to robert.danczak@pnnl.gov.

Here is a brief outline of the included files:

<b>Carbon_Efficiency.R:</b> This script will calculate Lamba and Y_csi, which help provide insight into how efficient microbes can use a given molecular formula. This will generate a single figure.<br>
<b>FTICR Analysis Walkthrough.docx:</b> This Word document has step-by-step instructions which will walk through how to analyze the provided FTICR data<br>
<b>FTMS_Analysis.R:</b> This script is a necessary pre-processing step for a provided FTICR report in order for other scripts to function. This will remove isotopic peaks and limit the mass range of peaks to 200-900 m/z. This script is based upon the <a href=https://github.com/EMSL-Computing/ftmsRanalysis/tree/master>'ftmsRanalysis' R package</a><br>
<b>HJ_Andrews_Tutorial_Report.csv:</b> An example FTICR report generated using Formularity from Watershed 1 within the HJ Andrews Experimental Forest, Oregon<br>
<b>SBR_SFA_ICR_Tutorial.pptx:</b> The PowerPoint shown at the Watershed Workshop ICR Tutorial 2019<br>
<b>Transformation_Analysis.R:</b> This script will perform a "transformation analysis" which will identify putative biochemical reactions within a dataset given a transformation database. 
Such analyses are outlined in <a href=https://doi.org/10.1007/s11306-006-0029-z>Breitling et al., 2006 - Metabolomics</a> and <a href=https://doi.org/10.1002/2017JG003967>Graham et al., 2017 - JGR: Biogeosciences</a>. This script will generate two figures.<br>
<b>Transformation_List_Amino_Acid.csv:</b> The transformation database for the transformation analysis; only identifies putative amino acid transformations.<br>
<b>install_packages_commands.txt:</b> A file to allow copying-and-pasting of R install package commands.<br>
<br>
The <b>AGU CZO:Watershed Workshop 2019</b> folder contains the necessary files for the FTICR tutorial given at that workshop.<br>
<br>
The <b> EMSL Summer School 2020 </b> folder contains all of the necessary files for the tutorial given during the summer school.
