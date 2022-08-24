# TFM_Alejandro_Santos
This repository aims to store the code and the data used for the analysis of the HDL proteome during postprandial hyperlipidemia. 

## Script "Proteomic_HDL_hyperlipidemia.R"
The Script can be seem as two parts, the first one focused in the preprocessing of the data. Then, the imputed and normalized dataset are used for the differential analysis by amica v2022.07.08 (check author repository https://github.com/tbaccata/amica) and statistical testing by MetaboAnalyst v5.0 (https://www.metaboanalyst.ca/). Amica results are then visualise in the Script. The second part of the Script focused on the coexpression clustering analysis using coseq v3.15 package (https://bioconductor.org/packages/release/bioc/html/coseq.html), in this part the dataset used is the raw data but with missing values imputed by the minimum value in each raw. Finally, the Script visualise the clustering results.

## Data 
In the Data folder it can be found all the needed files for the script to be run and also the results of the srcipt once it have been fully run. In the Amica folder, you will find the 4 required files to run the programa: contrast matrix, input, experimental desing and specification. Take into account that the input will me the same file as the final output of the preprocessing part of the Script (dietas_norm_imputed), this also applied for MetaboAnalyst input, but it was fragmented depending on the statistical test we want to perform. 
