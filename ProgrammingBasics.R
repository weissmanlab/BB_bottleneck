# Filename: ProgrammingBasics.R
 
# ---Simple Calculations---
2 + 3
 
x <- 2
y <- 3
x + y
x * y
 
# ---Data Structures---
 
# Vectors
workshop <- c(1, 2, 1, 2, 1, 2, 1, 2)
print(workshop)
workshop
 
gender <- c("f", "f", "f", NA, "m", "m", "m", "m")
q1 <- c(1, 2, 2, 3, 4, 5, 5, 4)
q2 <- c(1, 1, 2, 1, 5, 4, 3, 5)
q3 <- c(5, 4, 4,NA, 2, 5, 4, 5)
q4 <- c(1, 1, 3, 3, 4, 5, 4, 5)
 
# Selecting Elements of Vectors
q1[5]
q1[ c(5, 6, 7, 8) ]
q1[5:8]
q1[gender == "m"]
mean( q1[ gender == "m" ], na.rm = TRUE)
 
# ---Factors---
 
# Numeric Factors
 
# First, as a vector
workshop <- c(1, 2, 1, 2, 1, 2, 1, 2)
workshop
table(workshop)
mean(workshop)
gender[workshop == 2]
 
# Now as a factor
workshop <- c(1, 2, 1, 2, 1, 2, 1, 2)
workshop <- factor(workshop)
workshop
table(workshop)
mean(workshop) #generates error now.
gender[workshop == 2]
gender[workshop == "2"]
 
# Recreate workshop, making it a factor
# including levels that don't yet exist.
workshop <- c(1, 2, 1, 2, 1, 2, 1, 2)
workshop <- factor(
workshop,
levels = c( 1,   2,     3,      4),
labels = c("R", "SAS", "SPSS", "Stata")
)
 
# Recreate it with just the levels it
# curently has.
workshop <- c(1, 2, 1, 2, 1, 2, 1, 2)
workshop <- factor(
workshop,
levels = c( 1,  2),
labels = c("R","SAS")
)
 
workshop
table(workshop)
gender[workshop == 2]
gender[workshop == "2"]
gender[workshop == "SAS"]
 
# Character factors
 
gender <- c("f", "f", "f", NA, "m", "m", "m", "m")
gender <- factor(
gender,
levels = c("m",    "f"),
labels = c("Male", "Female")
)
 
gender
table(gender)
workshop[gender == "m"]
workshop[gender == "Male"]
 
# Recreate gender and make it a factor,
# keeping simpler m and f as labels.
gender <- c("f", "f", "f", NA, "m", "m", "m", "m")
gender <- factor(gender)
gender
 
# Data Frames
mydata <- data.frame(workshop, gender, q1, q2, q3, q4)
mydata
 
names(mydata)
row.names(mydata)
 
# Selecting components by index number
mydata[8, 6] #8th obs, 6th var
mydata[ , 6] #All obs, 6th var
mydata[ , 6][5:8] #6th var, obs 5:8
 
# Selecting components by name
mydata$q1
mydata$q1[5:8]
 
# Example renaming gender to sex while
# creating a data frame (left as a comment)
#
# mydata <- data.frame(workshop, sex = gender,
#   q1, q2, q3, q4)
 
# Matrices
 
# Creating from vectors
#mymatrix <- cbind(q1, q2, q3, q4)
#mymatrix
#dim(mymatrix)
 
# Creating from matrix function
# left as a comment so we keep
# version with names q1, q2...
#
# mymatrix <- matrix(
#   c(1, 1, 5, 1,
#     2, 1, 4, 1,
#     2, 2, 4, 3,
#     3, 1, NA,3,
#     4, 5, 2, 4,
#     5, 4, 5, 5,
#     5, 3, 4, 4,
#     4, 5, 5, 5),
#  nrow = 8, ncol = 4, byrow = TRUE)
# mymatrix
 
