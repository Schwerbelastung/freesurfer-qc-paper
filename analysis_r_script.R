#######################################################################################
### Script for investigating brain surface area properties between different groups ###
#######################################################################################

# Empty workspace
rm(list = ls())


#################################
### 1. Set up the environment ###
#################################

working_dir <- "/..."

setwd(working_dir) # Set working directory

# Load CSV files

unprocessed_csv <- paste(working_dir,"/thickness_BL_unprocessed.csv",sep="")
processed_csv <- paste(working_dir,"/thickness_BL_processed.csv",sep="")
additional_data <- paste(working_dir,"/additional_data.csv",sep="")
additional_data_gender <- paste(working_dir,"/additional_data_gender.csv",sep="")
additional_data_wilcoxon_sex <- paste(working_dir,"/additional_data_wilcoxon_sex.csv",sep="")
delimiter <- ";";

df_unprocessed <- read.csv(file = unprocessed_csv, sep=delimiter)
#df_unprocessed_original <- df_unprocessed # Store backup of original

df_processed <- read.csv(file = processed_csv, sep=delimiter)
#df_processed_original <- df_processed # Store backup of original

additional_data <- read.csv(file = additional_data, sep=delimiter)
#additional_data_original <- additional_data # Store backup of original

additional_data_gender <- read.csv(file = additional_data_gender, sep=delimiter)
#additional_data_original <- additional_data # Store backup of original

additional_data_wilcoxon_sex <- read.csv(file = additional_data_wilcoxon_sex, sep=delimiter)

df_bl_vs_fu_errors_grouped <- paste(working_dir,"/errors_bl_vs_fu_grouped.csv",sep="")
df_bl_vs_fu_errors_grouped <- read.csv(file = df_bl_vs_fu_errors_grouped, sep=delimiter)






###########################################
### 2. Calculate tables for differences ###
###########################################


df_differences_absolute = df_processed # Copy the table structure


for (i in 1:dim(df_processed)[1]) {
  for (j in 5:dim(df_processed)[2]) {
    value <- df_processed[i,j] - df_unprocessed[i,j]
    df_differences_absolute[i,j] = value
  }
}



df_differences_percentage = df_processed # Copy the table structure


for (i in 1:dim(df_processed)[1]) {
  for (j in 5:dim(df_processed)[2]) {
    value <- ((df_processed[i,j] - df_unprocessed[i,j]) / df_processed[i,j]) * 100
    df_differences_percentage[i,j] = value
  }
}





# Visualize different pairs' correlation

library("ggplot2")
library("ggpubr")


ggscatter(additional_data, x = "AGE", y = "Errors", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          cor.coef.name = "rho",
          xlab = "Age", ylab = "Number of errors")

t.test(additional_data["AGE"],additional_data["Errors"])


ggscatter(additional_data, x = "Patient", y = "Errors", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Patient (1 = yes)", ylab = "Number of errors")

t.test(additional_data["Patient"],additional_data["Errors"])

wilcox.test(Errors ~ Patient, data = additional_data, paired = FALSE)


ggscatter(additional_data, x = "SEX", y = "Errors", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "SEX (1 = male)", ylab = "Number of errors") +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15))

t.test(additional_data["SEX"],additional_data["Errors"])

####################################
### Correlation of error numbers ###
####################################

errornumbers_csv <- paste(working_dir,"/BL_and_FU_error_comparison.csv",sep="")
df_BLandFUerrors <- read.csv(file = errornumbers_csv, sep=delimiter)

ggscatter(df_BLandFUerrors, x = "BL.errors", y = "FU.errors", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "BL errors", ylab = "FU errors")

wilcox.test(errors ~ group, data = df_bl_vs_fu_errors_grouped, paired = TRUE)


#t.test(x = df_BLandFUerrors$BL.errors, 
#       y = df_BLandFUerrors$FU.errors, 
#       paired=TRUE)


############################################################
### Error number means and group analysis between groups ###
############################################################

all_errors_csv <- paste(working_dir,"/all_errors_per_subject.csv",sep="")
df_errorspersubject <-  read.csv(file = all_errors_csv, sep=delimiter)

wilcox.test(Errors ~ Group, data = df_errorspersubject, paired = FALSE)

wilcox.test(errors ~ sex, data = additional_data_wilcoxon_sex, paired = FALSE)

mean(na.omit(df_errorspersubject$Errors_FEP))
sd(na.omit(df_errorspersubject$Errors_FEP))

mean(na.omit(df_errorspersubject$Errors_CTRL))
sd(na.omit(df_errorspersubject$Errors_CTRL))

################
### Barplots ###
################

# library
library(ggplot2)

# create a dataset
NumberOfErrors <- c(rep("0" , 2) , rep("01" , 2) , rep("02" , 2) , rep("03" , 2),
            rep("04" , 2) , rep("05" , 2) , rep("06" , 2) , rep("07" , 2),
            rep("08" , 2) , rep("09" , 2) , rep("10" , 2) , rep("11" , 2),
            rep("12" , 2) , rep("13" , 2) , rep("14" , 2) )
Sex <- rep(c("male" , "female") , 15)
Errors <- c(2,11,13,10,8,4,10,3,12,2,5,2,9,2,2,0,3,0,1,0,0,0,1,0,1,0,1,0,0,0)#abs(rnorm(12 , 0 , 15)) #double-checked
data <- data.frame(NumberOfErrors,Sex,Errors)

# Grouped
ggplot(data, aes(fill=Sex, y=Errors, x=NumberOfErrors)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  #ggtitle("Image error counts") +
  xlab("Error count in individual images") +
  ylab("Number of images in dataset") +
  theme(plot.title = element_text(hjust = 0.5))
  


########