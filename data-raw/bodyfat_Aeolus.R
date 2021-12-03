## code to prepare `bodyfat_Aeolus` dataset goes here

data <- read.csv("data-raw/ChengJAE2018_data4dryad.csv", header = TRUE, sep = ";", dec = ".")
data <- data[, -c(1,4,6,7,8,9,10)]

year  <- c(rep("2016", 239), rep("2009", 350))
data  <- cbind(data, year)

data <- na.omit(data)
rownames(data) <- NULL

aux <- data[data$site == "Aeolus", -c(1)]
rownames(aux) <- NULL

days <- NULL
for(i in 1:length(aux[,1])){
  if( aux$date[i] == "2015-Nov-17"){
    days[i] <- 127
  }
  if( aux$date[i] == "2015-Nov-13"){
    days[i] <- 130
  }
  if( aux$date[i] == "2016-Mar-09"){
    days[i] <- 15
  }
  if( aux$date[i] == "2008-Nov-18"){
    days[i] <- 126
  }
  if( aux$date[i] == "2009-Jan-31"){
    days[i] <- 52
  }
  if( aux$date[i] == "2009-Mar-27"){
    days[i] <- 0
  }
  if( aux$date[i] == "2015-Nov-05"){
    days[i] <- 139
  }
  if( aux$date[i] == "2016-Mar-10"){
    days[i] <- 13
  }
  if( aux$date[i] == "2008-Nov-18"){
    days[i] <- 126
  }
  if( aux$date[i] == "2009-Jan-12"){
    days[i] <- 71
  }
  if( aux$date[i] == "2009-Mar-19"){
    days[i] <- 4
  }
  if( aux$date[i] == "2015-Nov-06"){
    days[i] <- 138
  }
  if( aux$date[i] == "2015-Nov-05"){
    days[i] <- 139
  }
  if( aux$date[i] == "2016-Mar-04"){
    days[i] <- 20
  }
  if( aux$date[i] == "2016-Feb-17"){
    days[i] <- 35
  }
  if( aux$date[i] == "2008-Feb-12"){
    days[i] <- 40
  }
  if( aux$date[i] == "2008-Mar-26"){
    days[i] <- 0
  }
  if( aux$date[i] == "2015-Nov-18"){
    days[i] <- 126
  }
  if( aux$date[i] == "2016-Mar-21"){
    days[i] <- 2
  }
  if( aux$date[i] == "2008-Nov-25"){
    days[i] <- 119
  }
  if( aux$date[i] == "2009-Jan-13"){
    days[i] <- 70
  }
  if( aux$date[i] == "2009-Mar-15"){
    days[i] <- 8
  }
  if( aux$date[i] == "2009-Mar-16"){
    days[i] <- 7
  }
  if( aux$date[i] == "2008-Dec-17"){
    days[i] <- 97
  }
}

days <- 175 - days
bodyfat_Aeolus <- cbind(aux[,-1], days)
colnames(bodyfat_Aeolus) <- c("sex", "percentfat", "year", "days")

bodyfat_Aeolus$sex <- factor(bodyfat_Aeolus$sex)
bodyfat_Aeolus$year <- factor(bodyfat_Aeolus$year)

usethis::use_data(bodyfat_Aeolus, overwrite = TRUE)
