#############################################################################################
#### Visualizing template distribution and setting the minimum and maximum window length ####
#############################################################################################

setwd("/local/path/")

#read the txt file that was produced in 05.1. FragmentLengthRange.slurm
template.length=read.table("C:/Users/viode560/Documents/statistics_work/template_size/template_length_HPT.txt", sep="\t", header=FALSE)

# visualization in the form of an histogram
hist(abs(template.length$V1), main = "Template length hypothalamus", xlab = "Base pairs")

#Assuming template.length is a data frame, the template.length$V1 syntax extracts the V1 column from the template.length data frame. 
#The abs() function takes the absolute value of each element in the resulting column. 
#The hist() function then creates a histogram of the absolute values of V1, with the specified title and x-axis label.

#summary to have the minimum and maximum window length
stats.HPT=summary(abs(template.length$V1))
stats.HPT
