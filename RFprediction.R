.libPaths(c('R',.libPaths()))

#require(randomForest) || 
#install.packages("randomForest", repos="http://cran.r-project.org")
library(randomForest,lib.loc = "/WWW/httpd-users/molgen/Sorek/PASIFIC/R")

# load the random forest 
load("rf13-3.RData")

# get out filename
file <- commandArgs(trailingOnly = TRUE)[1]

# load user's data
userdata = read.delim(file, header = FALSE, sep="\t", col.names = c("ID", "seq", 
                                                                      "TM", "TM.seq", "start34", "len34", "end34", "G34", "G34.len", "Gfold34", "Gfold34.len", "freq_MFE34", "diversity34", "avg_prob34", "avg_sl34_prob", "mul_prob34", "mul_sl34_prob", 
                                                                      "AT", "AT.seq", "start23", "len23", "end23", "G23", "G23.len", "Gfold23", "Gfold23.len", "freq_MFE23", "diversity23", "avg_prob23", "avg_sl23_prob", "mul_prob23", "mul_sl23_prob", 
                                                                      "P1", "P1.seq", "start", "len", "end", "GP1", "GP1.len", 
                                                                      "overlap", "ddG", "fcG", "ddG.SL", "fcG.SL", 
                                                                      "G", "fr", "di", 
                                                                      "Gfold23.G", "Gfold23.G.1", "Gfold34.G", "Gfold34.G.1", 
                                                                      "dprob", "dsl_prob", "d_mul_prob", "dsl_mul_prob",
                                                                      "TMfold", "ATfold"))


# Predict
Prediction <- predict(fit, userdata, type="prob")

notes <- userdata$TM
notes <- sub("[.()]+", "", notes)

# write to CSV file
submit = data.frame(ID = userdata$ID, 
                    Score = Prediction[,2], 
                    Seq = userdata$seq, 
                    Notes = notes,
                    TMfold = userdata$TMfold, 
                    TMfold_free_energy = userdata$Gfold34,
                    ATfold = userdata$ATfold,
                    ATfold_free_energy = userdata$Gfold23)
write.csv(submit, file = paste(file, ".prediction.csv", sep=""), row.names = FALSE)