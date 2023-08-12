
#-------------------------|-----------------------------------------------------
# Packages & Dependencies | 
#-------------------------|-----------------------------------------------------


library(readxl)       # for reading xlxs
library(stringr)      # for parsing strings
library(dplyr)        # for data frame manipulation
library(latticeExtra) # for visualizations
library(ggplot2)      # for visualizations
library(treemap)      # for treemaps
library(xtable)       # for dataframe to latex

library(plotly)       # for 3D PCA Plot
library(MASS)         # for LDA/QDA
library(fda)          # for FLDA/FQDA
library(klaR)         # for FLDA


library(corrplot)     # for correlation plots
library(factoextra)   # for PCA visualization

library(e1071)        # for svm analysis
library(caret)        # for confusionMatrix()

#------------------------------------------------------------------------------
#               Import Data  
#------------------------------------------------------------------------------

DD<-read_xlsx("ssc_case_study_2_inflammatory_bowel_disease.xlsx", sheet=1)
DD<-as.data.frame(DD)

#------------------------------------------------------------------------------
#               Data Wrangling 
#------------------------------------------------------------------------------

DF<-DD
DF[1:5,1:6]

# ---Change Column Names ---

# we want to shorten each subject ID and remove redundant information
class(colnames(DF))
ID<-colnames(DF)

#reduce ID to be the unique ID of each individual
id<-str_extract(ID[3:length(ID)], "(?<=_)\\d+")

#make sure each ID is unique
is_unique <- sum(duplicated(id)) == 0
print(is_unique)

#assigning the new `id` to the column names 
colnames(DF)[3:length(ID)]<-id

dim(DF) # 313 rows, 128 cols


#------CHECK FOR NA VALUES------
sum(is.na(DF)) # 313 values in total
which(is.na(DF))
#NA values of the first column
sum(is.na(DF[,1])) # 4 NA Values in col 1
which(is.na(DF[,1])) # first 4 rows of col 1 is NA
#NA values of the second column
sum(is.na(DF[,2])) #309 NA Values in col 2
which(!is.na(DF[,2])) # only first 4 rows of col 2 are not NAs

# rest of the columns
sum(is.na(DF[,3:128])) # no other NA values

# we combine the first 2 columns to eliminate NA values

var<-c()
var[1:4]<-DF[1:4,2]
var[5:313]<-DF[5:313,1]
head(var)
V<-as.data.frame(var)

df<-cbind(V,DF[,3:128])

#preview dataframe first 10 column, 10 rows
df[1:10,1:10]

#------TRANSPOSE DATAFRAME------

# we combine `probe set id` and other variables columns into one
# by doing so we eliminate the NA values

tdf<-t(as.matrix(df))
df<-as.data.frame(tdf)
df[1:5,1:5]              # preview data
dim(df)                  # new dimensions are 127 rows 313 cols

names<-c(df[1,])         # save first row as names of cols
colnames(df)<-names      # change col names

dim(df)

df<-df[2:127,]           # remove first row

# subjectID<-row.names(df) # save subject ID
# df<-cbind(subjectID,df)  # adding subjectID as column variable 

dim(df) # our final dimensions after manipulation are 126x313

class(df$`200006_at`) # we notice that the Probe Set Values are `character`
all(sapply(df[,5:313], is.character))


#------WRANGLE VARIABLES------

# Ethnicity -> As Factor
unique(df$Ethnicity)
# correcting mis-spelling
sum(df$Ethnicity=="cacuasian") # caucasion mis-spelled once
which(df$Ethnicity=="cacuasian") # index 97

df$Ethnicity[97]<-"caucasian"

df$Ethnicity<-as.factor(df$Ethnicity)

#------Convert Age to Numeric------
class(df$Age)
df$Age<-as.numeric(df$Age)

# create age group variable

max(df$Age)      # 73 Years 
min(df$Age)      # 20 Years

# Creating Age Class

# 20 =< age < 35    "A: Adult"              df$Age<35
# 35 =< age < 50    "B: Middle Age"         df$Age>=35 & df$Age<50
# 50 =< age         "C: Senior"             df$Age>=50


df<-df%>%mutate(AgeClass = ifelse(Age<35,"A",
                              ifelse(Age>=35 & Age<50, "B", "C")))

#------Convert Sex to Binary------

# 1 -> MALE
# 0 -> FEMALE

df$Sex<-recode(df$Sex, "female"=0, "male"=1)

#------Convert Group (Disease)------

# correcting mis-spelling
unique(df$Group) 
sum(df$Group=="Ulcerative") # 1 observation
which(df$Group=="Ulcerative") # index 91st observation
df$Group[91]<-"Ulcerative Colitis"

# Coding Disease into 3 Categories
# Levels:  NN -> Normal, UC -> Ulcerative Colitis,  CD -> Crohn's Disease


df$Group<-recode(df$Group, "Ulcerative Colitis"="UC", 
                            "Crohn's Disease"="CD", 
                            "Normal"="NN")

df$Group<-as.factor(df$Group) # changing diseases as factors
df[1:8,1:8]

df<-cbind((df[,1:4]),(df[,314]),(df[,5:313]))
colnames(df)[5]<-"AgeClass"
#------Convert Probe Set Values to Numeric------

df<- df%>% mutate_at(vars(6:314), as.numeric) # converts columns 6: to numeric

#------Save Wrangled Dataset-----

# write.csv(df,file = "Gene.csv", row.names = TRUE)



#------------------------------------------------------------------------------
#               EXPLORATORY DATA ANALYSIS
#------------------------------------------------------------------------------


summary(df[,2:6])
table(df$Group)
table(df$Sex)
table(df$Ethnicity)
summary(df)

table(df$Sex,df$Group)

# Histogram of number of subjects in study by ethnicity
# Possible Visualizations in Tableau

df<-read.csv("Gene.csv")


#------Age Distribution by Grouped Age Class ------

# Histogram by age groups
ggplot(df, aes(x = Age, fill = cut(Age, 
                                   breaks = c(19, 35, 50, 65, Inf), 
                                   labels = c("20-35", "35-50", "50-65", "65+")
                                   )
)) + geom_histogram(binwidth = 5, position = "dodge", 
                    color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("#E2B6C7", "#B8A9C9", "#98B6B1", "#B6C9A8")) +
  labs(x = "Age", y = "Frequency", fill = "Age Class") +
  theme_classic()


# Density of Age Distribution fill by age groups
ggplot(df, aes(x = Age,fill="#98B6B1")) + 
  geom_density(color = "black", alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = "#98B6B1") +
  labs(x = "Age", y = "Proportion") + ggtitle("Age Distribution")+
  theme_classic()

# Density of Age Distribution fill by age groups
ggplot(df, aes(x = Age, fill = cut(Age, 
                                   breaks = c(19, 35, 50, 65, Inf), 
                                   labels = c("20-35", "35-50", "50-65", "65+")
)
)) + geom_density(color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("#E2B6C7", "#B8A9C9", "#98B6B1", "#B6C9A8")) +
  labs(x = "Age", y = "Proportion", fill = "Age Class") + 
  ggtitle("Age Distribution by Age Class")+
  theme_classic()


#------Ethnicity ------
# Ethinicity Tree Map

eth_freq<-table(df$Ethnicity)
ggplot(df, aes(x = Ethnicity)) +
  geom_bar(fill = "#98B6B1") +
  labs(x = "Ethnicity", y = "Frequency") +
  ggtitle("Frequency of Ethnicity") +
  theme_classic()
ggplot(df, aes(x = Ethnicity, fill = Ethnicity)) +
  geom_bar() +
  scale_fill_manual(values =  c("#B39EB5", "#8FBEB7", 
                                "#B4D4AB", "#E4C1B7", "#F3A0B3")) +
  labs(x = "Ethnicity", y = "Frequency") +
  ggtitle("Frequency of Ethnicity") +
  theme_classic()


# Create a data frame with the proportion of each level of "group"
group_prop <- df %>% 
  group_by(Group) %>% 
  summarize(prop = n()/nrow(df))

group_prop$percent <- scales::percent(group_prop$prop)
# Disease Proportion
ggplot(group_prop, aes(x = "", y = prop, fill = Group)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9DDCDC", "#C4C4C4","#FFB6C1")) +
  labs(fill = "Group")+ ggtitle("Disease Status")+
    geom_text(aes(label = percent), position = position_stack(vjust = 0.5))




#------ Sex ------

sex_count <- table(df$Sex)
df_sex_count <- data.frame(Sex = names(sex_count), 
                           Count = as.numeric(sex_count))

ggplot(df_sex_count, aes(x = Sex, y = Count, fill = Sex)) +
  geom_col() +
  scale_fill_manual(values = c("#98B6B1", "#F6AE2D"), 
                    labels=c("Female","Male")) +
  labs(x = "Sex", y = "Count", fill = "") + ggtitle("Sex Distribution")+
  theme_minimal()
#------------------------------------------------------------------------------
#               DIMENSION REDUCTION
#------------------------------------------------------------------------------

#------ PCA Principle Component Analysis ------
# we will use principle compenent analysis to reduce the data into 
# 3 principle axes, this way we can visualize the high dimensional
# data in a better way



#PRINCIPLE COMPONENT ANALYSIS

# # Variance Matrix
# nn<-df[,5:313] # subset dataframe to only numeric columns


nn<-df[,-c(3,5)]
S<-var(nn[2:312])     # compute variance matrix S

 
# Eigen Values & Eigen Vectors

eigenS<-eigen(S)
Eval<-eigenS$values
Evec<-eigenS$vectors
# adjoint(eigenS$vectors)
# print(eigenS, digits=5)


# #Proportion of total variance due to (explained by) the kth principal 
# component is
# # First Principle Component
p1<-eigenS$values[1]/sum(diag(S))
# # First Two Principle Components
p2<-(eigenS$values[1]+eigenS$values[2])/sum(diag(S))
# 
# # First Three Principle Components
p3<-(eigenS$values[1]+eigenS$values[2]+eigenS$values[3])/sum(diag(S))




# Extract the numeric data from your data frame

pcadata <- nn[,-c(1)]
pcadata <- scale(pcadata)  # Scale the numeric data

# Perform PCA
pca <- prcomp(pcadata)

summary(pca)


# Extract the first three principal components
prc <- pca$x[, 1:3]

# Add the group labels back to the data
datapca <- cbind(prc, Group = nn$Group)
datapca<-as.data.frame(datapca)
# Define a color palette for the groups
colors<-c("#FFB6C1","#87CEFA","#90EE90")
# Create a 3D scatter plot with Plotly
p <- plot_ly(datapca, x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~Group, colors = colors, 
             type = "scatter3d", mode = "markers")

# Add axis labels and a title to the plot
p <- layout(p, scene = list(xaxis = list(title = "PC1"), 
                            yaxis = list(title = "PC2"), 
                            zaxis = list(title = "PC3")), 
            title = "PCA Plot")

# Show the plot
p


# proportions of variation explained

proportions<-as.matrix(pca$sdev^2/sum(pca$sdev^2))
cumulative<-cumsum(proportions)

cbind(proportions,cumulative)

# looking at the first 3 components

cumulative[3] # captures 0.5584819% of variation in the data

# which is not very good.

# to obtain
cumulative[8]  # 0.8135988 % captured by first 8 principle components
cumulative[14] # 0.9072991 % captured by first 14 principle components

#PCA reduces our dimensionality down to 14

# #------ LDA Fisher's Linear Discriminant Analysis ------


#------------------------------------------------------------------------------
#               Binary Classification
#------------------------------------------------------------------------------

# looking at the binary classification problem, we focus on UC and CD patients
# this examines the difference between the two types of diseases

# UC vs CD
bf <- subset(df, Group %in% c("UC", "CD"))
bf$Group<- as.character(bf$Group)
bf$Group<- as.factor(bf$Group)
levels(bf$Group)

bf$disease <- ifelse(bf$Group == "UC", 1, 0)

bf<-bf[,2:ncol(bf)]

# Dimension Reduction to 1 dimension





ppCD<-as.matrix(table(bf$disease))[1]/nrow(bf) #prior of CD
ppUC<-as.matrix(table(bf$disease))[2]/nrow(bf) #prior of UC

priors<-c(ppCD,ppUC)
priors

#------------------------------------------------------------------------------
#               CROSS - VALIDATION
#------------------------------------------------------------------------------

# We have 126 Patients
# We will use a 80/20 Data Split

#------------------------------------------------------------------------------
#               CLASSIFICATION -- SVM W/O  Dimension Reduction
#------------------------------------------------------------------------------

# Test/Train Split

#dim(df) : 126x314
df<-df[,-c(5)]# leaves out column 5 (age class)
# we remove the age class column as we only created it for EDA purposes.
# it is created from `Age`, thus including both age and age class in the model 
# would be redundant


# df is now a dataframe of 313 columns, after exluding age class
nn<-df[,-c(3)] # only numerical predictors (excludes ethnicity)

dim(nn) # nn is a dataframe of 312 columns after exluding ethinicity

set.seed(1)
# index<-sample(seq(1,126,1),size=84,replace = FALSE)
index<-sample(seq(1,126,1),size=round(126*0.8),replace = FALSE)

train<-df[index,]
test<-df[-index,]

dim(test)
dim(train)


svm_model = svm(Group~., data = train, 
                kernel = "poly",
                type="C-classification",
                gamma=0.1,
                cost = 10, 
                scale = FALSE)
# summary(svm_model)
# length(svm_pred)
svm_pred<-predict(svm_model,newdata=test)
table(svm_pred,test$Group)
acc<-sum(diag(table(svm_pred,test$Group)))/25
confusionMatrix(data=svm_pred,test$Group)


#_______________________________________________

# SVM With K-Fold CV

set.seed(1)
n<-126
nt<-round(0.8*n) # 80/20 Split
neval<-n-nt
rep<-1000

errsvm<-dim(rep)
accsvm<-dim(rep)
for( k in 1:rep){
  id<-sample(1:n,nt, replace = FALSE)
  train<-df[id,]
  test<-df[-id,]
  
  #SVM
  svm_model = svm(Group~., data = train, 
                  kernel = "poly",
                  type="C-classification",
                  gamma=0.1, # low gamma vs high gamma
                  cost = 1,  # low cost vs high cost
                  scale = FALSE)
  svm_pred<-predict(svm_model,newdata=test)
  tab_svm<-table(svm_pred,test$Group)
  acc_svm<-sum(diag(table(svm_pred,test$Group)))/neval
  CM<-confusionMatrix(data=svm_pred,test$Group)
  errsvm[k]<-(neval-sum(diag(tab_svm)))/neval
  accsvm[k]<-acc_svm
}

MEsvm<-mean(errsvm) # mean of error_svm
MAsvm<-mean(accsvm) # mean of accuracy_svm
MEsvm 
MAsvm

summary(svm_model)
# Poly Kernel
# > MEsvm
# [1] 0.17468
# > MAsvm
# [1] 0.82532

# Linear Kernel
# > MEsvm
# [1] 0.18484
# > MAsvm
# [1] 0.81516

# Sigmoid Kernel
# > MEsvm
# [1] 0.5356
# > MAsvm
# [1] 0.4644

# Radial Kernel
# > MEsvm
# [1] 0.5356
# > MAsvm
# [1] 0.4644


#------ SVM after LDA ------


#Prior Probablities of Class NN, UC, and CD
table(df$Group) # table of frequencies for each class of Disease

as.matrix(table(df$Group))[1] # freq of CD
as.matrix(table(df$Group))[2] # freq of NN
as.matrix(table(df$Group))[3] # freq of UC

ppCD<-as.matrix(table(df$Group))[1]/nrow(df) #prior of CD
ppNN<-as.matrix(table(df$Group))[2]/nrow(df) #prior of NN
ppUC<-as.matrix(table(df$Group))[3]/nrow(df) #prior of UC
priors<-c(ppCD,ppNN,ppUC)
sum(priors)                 # make sure sum of priors are equal to 1
priors

LDT<-lda(Group~., data=nn, prior=priors)
dfLDA<-predict(LDT, nn)$x

dd<-as.data.frame(cbind(as.character(df$Group),dfLDA))
colnames(dd)<-c("Group","LD1","LD2")

dd$LD1<-as.numeric(dd$LD1)
dd$LD2<-as.numeric(dd$LD2)
dd$Group<-as.factor(dd$Group)


# Plot LDA results
P<-ggplot(data = dd, aes(x = LD1, y = LD2, color = Group)) +
  geom_point() +
  ggtitle("LDA Plot") +
  xlab("LD1") +
  ylab("LD2") +
  theme_bw()
P

# We observe from the plot that LD1 has good between 
# seperation between classes CD and NN
# the LD2 Dimension has good between seperation between CD and UC 
# the LD2 Dimension has good between seperation between NN and UC 

########


# LDA -> SVM With K-Fold CV Polynomial Kernal

set.seed(1)

n<-126
nt<-round(0.8*n) # 80/20 Split
neval<-n-nt
rep<-1000

errsvmx<-dim(rep)
accsvmx<-dim(rep)


for( k in 1:rep){
  idx<-sample(1:n,nt, replace = FALSE)
  trainx<-dd[idx, ]
  testx<-dd[-idx, ]
  
  #SVM
  svm_modelx = svm(trainx$Group~., data = trainx, 
                   kernel = "poly",
                   type="C-classification",
                   gamma=0.1, # low gamma vs high gamma
                   cost = 1,  # low cost vs high cost
                   scale = FALSE)
  svm_predx<-predict(svm_modelx,newdata=testx)
  tab_svmx<-table(svm_predx,testx$Group)
  acc_svmx<-sum(diag(table(svm_predx,testx$Group)))/neval
  CMx<-confusionMatrix(data=svm_predx,testx$Group)
  errsvmx[k]<-(neval-sum(diag(tab_svmx)))/neval
  accsvmx[k]<-acc_svmx
}

MEsvmx<-mean(errsvmx) # mean of error_svm
MAsvmx<-mean(accsvmx) # mean of accuracy_svm
MEsvmx #0.09408
MAsvmx #0.90592


########

# LDA -> SVM With K-Fold CV Linear Kernal

set.seed(1)
n<-126
nt<-round(0.8*n) # 80/20 Split
neval<-n-nt
rep<-1000

errsvmx<-dim(rep)
accsvmx<-dim(rep)
confusion<- list()

for( k in 1:rep){
  idx<-sample(1:n,nt, replace = FALSE)
  trainx<-dd[idx, ]
  testx<-dd[-idx, ]
  
  #SVM
  svm_modelx = svm(trainx$Group~., data = trainx, 
                   kernel = "linear",
                   type="C-classification",
                   gamma=0.1, # low gamma vs high gamma
                   cost = 1,  # low cost vs high cost
                   scale = FALSE)
  svm_predx<-predict(svm_modelx,newdata=testx)
  tab_svmx<-table(svm_predx,testx$Group)
  acc_svmx<-sum(diag(table(svm_predx,testx$Group)))/neval
  CMx<-confusionMatrix(data=svm_predx,testx$Group)
  
  errsvmx[k]<-(neval-sum(diag(tab_svmx)))/neval
  accsvmx[k]<-acc_svmx
}

CMx
MEsvmx<-mean(errsvmx) # mean of error_svm
MAsvmx<-mean(accsvmx) # mean of accuracy_svm
MEsvmx #0.0308
MAsvmx #0.9692


# Plot LDA SVM

# plot the SVM model
plot(svm_modelx, trainx)




# summary table

UCs<-df %>% filter(Group == "UC")%>%
  summarize( samples = sum(Group=="UC"),
             agemean = mean(age),
             agemed = median(age),
             num_male = sum(Sex == 1),
             num_female = sum(Sex == 0))


CDs<-df %>% filter(Group == "CD")%>%
  summarize( samples = sum(Group=="CD"),
             agemean = mean(age),
             agemed = median(age),
             num_male = sum(Sex == 1),
             num_female = sum(Sex == 0))

NNs<-df %>% filter(Group == "NN")%>%
  summarize( samples = sum(Group=="NN"),
             agemean = mean(age),
             agemed = median(age),
             num_male = sum(Sex == 1),
             num_female = sum(Sex == 0))


total<-df %>% summarize( samples = nrow(df),
             agemean = mean(age),
             agemed = median(age),
             num_male = sum(Sex == 1),
             num_female = sum(Sex == 0))
total

dtab<-cbind(t(CDs),t(UCs),t(NNs), t(total))
colnames(dtab)<-c("CD","UC","NN", "Total")
dtab
dtab<-as.data.frame(dtab)

print(xtable(dtab, type = "latex"), file = "ztex.tex")




