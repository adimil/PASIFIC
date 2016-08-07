# PASIFIC
http://www.weizmann.ac.il/molgen/Sorek/PASIFIC/

This is the source code for **PASIFIC - Prediction of Alternative Structures for Identification of *Cis*-regulation**.

##The pipeline:
#### 1. Run **structure_prediction.pl** for prediction of the alternative structures:

###### perl structure_prediction.pl \<tmpdir\> \<run_name\> \<p_P1\> \<p_simpAT\> \<p_limitAT\> \<p_findterm\>
* \<tmpdir\> - The directory for tem files (fasta should also be there)
* \<run_name\> - The name of the FASTA file (without .fasta extention)
* \<p_P1\> - Look for P1? (1/0)
* \<p_simpAT\> - Should the terminator be simple? (1/0)
* \<p_limitAT\> - Limit antiterminator length? (1/0)
* \<p_findterm\> - Are the TSSs unknown? (1/0)

#### 2. Run **RFprediction.R** for classification scores:
###### Rscript RFprediction.R \<struct_file.out\>
* \<struct_file.out\> - The perl script output file
* **rf13-3.RData** is the built classifier. Should reside in the same folder.

## Dependencies:
- perl
- The ViennaRNA Package (https://www.tbi.univie.ac.at/RNA/)
- Rscript
- randomForest R-package (https://cran.r-project.org/web/packages/randomForest/)

##Build classifier
Use **BuildClassifier.R** to rebuild the classifier.
Data available in:
**alldata.train.csv** - Data that were used for training.
**alldata.test.csv** - Data that were used for testing.
