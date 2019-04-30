#Data analysis for the Ising model
##################################

library(ggplot2)

myPalette<-c("black", "red")

#Some useful functions
######################
readData<-function(algorithmChoice){
  fileNames <- list.files(path="./Data", pattern=paste("A=", algorithmChoice, sep=""))
  dataSet <- lapply(paste("./Data/", fileNames, sep=""), read.csv)
  df <- do.call(rbind, dataSet)
  return(df)
}

readData3D<-function(){
  fileNames <- list.files(path="./SARW3D/Data", pattern="3D")
  dataSet <- lapply(paste("./SARW3D/Data/", fileNames, sep=""), read.csv)
  df <- do.call(rbind, dataSet)
  return(df)
}

readThermalizationTimes<-function(algorithmChoice, method){
  fileNames <- list.files(path="./Data/Thermalization Times/", 
                          pattern=paste("part",method,".*A=", algorithmChoice, sep=""))  
  dataSet <- sapply(paste("./Data/Thermalization Times/", fileNames, sep=""), read.table)
  usedSARWThermalizationLengths<-50*c(1:5)+1
  dataNames<-names(dataSet)
  for(i in 1:length(dataSet)){
    names(dataSet)[grep(paste("S=",usedSARWThermalizationLengths[i], sep=""), 
                        dataNames, value=FALSE)]<-paste("L", usedSARWThermalizationLengths[i], sep="")
  }
  return(dataSet)
}

readThermalizationTimes<-function(){
  dataSet<-list()
  dataSet[[1]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=51_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[2]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=101_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[3]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=151_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[4]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=201_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[5]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=251_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[6]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=51_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[7]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=101_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[8]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=151_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[9]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=201_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[10]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=251_A=Snake.txt", quote="\"", header=FALSE)
  dataSet[[11]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=51_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[12]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=101_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[13]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=151_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[14]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=201_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[15]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=251_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[16]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=51_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[17]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=101_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[18]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=151_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[19]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=201_A=Pivot.txt", quote="\"", header=FALSE)
  dataSet[[20]]<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=251_A=Pivot.txt", quote="\"", header=FALSE)
  
  
####  
  thermalizationTimeA51Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=51_A=Snake.txt", quote="\"")
  thermalizationTimeA101Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=101_A=Snake.txt", quote="\"")
  thermalizationTimeA151Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=151_A=Snake.txt", quote="\"")
  thermalizationTimeA201Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=201_A=Snake.txt", quote="\"")
  thermalizationTimeA251Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=251_A=Snake.txt", quote="\"")
  thermalizationTimeB51Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=51_A=Snake.txt", quote="\"")
  thermalizationTimeB101Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=101_A=Snake.txt", quote="\"")
  thermalizationTimeB151Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=151_A=Snake.txt", quote="\"")
  thermalizationTimeB201Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=201_A=Snake.txt", quote="\"")
  thermalizationTimeB251Snake<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=251_A=Snake.txt", quote="\"")
  thermalizationTimeA51Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=51_A=Pivot.txt", quote="\"")
  thermalizationTimeA101Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=101_A=Pivot.txt", quote="\"")
  thermalizationTimeA151Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=151_A=Pivot.txt", quote="\"")
  thermalizationTimeA201Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=201_A=Pivot.txt", quote="\"")
  thermalizationTimeA251Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partA_S=251_A=Pivot.txt", quote="\"")
  thermalizationTimeB51Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=51_A=Pivot.txt", quote="\"")
  thermalizationTimeB101Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=101_A=Pivot.txt", quote="\"")
  thermalizationTimeB151Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=151_A=Pivot.txt", quote="\"")
  thermalizationTimeB201Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=201_A=Pivot.txt", quote="\"")
  thermalizationTimeB251Pivot<-read.table("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/Thermalization Times/Thermalization_Time_partB_S=251_A=Pivot.txt", quote="\"")
  return(dataSet)
}

computeThermalizationTime<-function(data){
  thermalizationTime<-seq(0, 0, length=length(unique(data$Iteration)))
  
  for(i in 1:length(data$Iteration)){
    thermalizationTime[data$Iteration[i]]<-thermalizationTime[data$Iteration[i]]+1
  }
  return(thermalizationTime)
}

computeAverageTermalizationTime<-function(data){
  averageTermalizationTime<-{}
  
  return()
}

computeEndToEndAutocorrelation<-function(data){
  nSamples<-length(data$End_to_end_length)
  autocorrelation<-{}
  for(i in 1:(nSamples/2)){
    autocorrelation[i]<-sum((data$End_to_end_length[1:(nSamples/2)]-mean(data$End_to_end_length))*
                              (data$End_to_end_length[i:((nSamples/2)+i-1)]-mean(data$End_to_end_length)))
  }
  autocorrelation<-autocorrelation/autocorrelation[1]
  return(autocorrelation)
}

computeAverageEndToEndLength<-function(data, usedLengths){
  averageEndToEndLength<-{}
  for(i in 1:length(usedLengths)){
    averageEndToEndLength[i]<-mean(data[data$Size==usedLengths[i],]$End_to_end_length)
  }
  return(averageEndToEndLength)  
}

computeVarianceEndToEndLength<-function(data, usedLengths){
  varianceEndToEndLength<-{}
  for(i in 1:length(usedLengths)){
    varianceEndToEndLength[i]<-var(data[data$Size==usedLengths[i],]$End_to_end_length)
  }
  return(varianceEndToEndLength)  
}

computeAverageGyrationRadius<-function(data, usedLengths){
  averageGyrationRadius<-{}
  for(i in 1:length(usedLengths)){
    averageGyrationRadius[i]<-mean(data[data$Size==usedLengths[i],]$Gyration_radius)
  }
  return(averageGyrationRadius)  
}

computeExponent<-function(averageValues, usedLengths){
  linearFit<-lm(log(averageValues)~log(usedLengths))
  return(as.vector(linearFit$coefficients))
}

computeCorrelationTime<-function(autocorrelation, linearPartLength=Inf){
  endPointCorrelation<-min(which(autocorrelation<=0)-1)
  endPointLinearPartCorrelation<-min(endPointCorrelation, linearPartLength)
  autocorrelation.df<-data.frame(y=autocorrelation[1:endPointCorrelation],
                                              x=seq(endPointCorrelation))
  
  linearFitCoefficients<-as.vector(coef(lm(log(autocorrelation.df[,1])~autocorrelation.df[,2], 
                                           subset=c(1:endPointLinearPartCorrelation))))
  
  correlationTime<--1/linearFitCoefficients[2]
  return(correlationTime)
}
####


#Start the clock!
time<-proc.time()

# #Slithering Snake thermalization data
# #####################################
# thermalizationDataA51Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=51_A=2.txt")
# thermalizationTimeA51Snake<-computeThermalizationTime(thermalizationDataA51Snake)
# rm(thermalizationDataA51Snake)
# thermalizationDataA101Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=101_A=2.txt")
# thermalizationTimeA101Snake<-computeThermalizationTime(thermalizationDataA101Snake)
# rm(thermalizationDataA101Snake)
# thermalizationDataA151Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=151_A=2.txt")
# thermalizationTimeA151Snake<-computeThermalizationTime(thermalizationDataA151Snake)
# rm(thermalizationDataA151Snake)
# thermalizationDataA201Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=201_A=2.txt")
# thermalizationTimeA201Snake<-computeThermalizationTime(thermalizationDataA201Snake)
# rm(thermalizationDataA201Snake)
# thermalizationDataA251Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=251_A=2.txt")
# thermalizationTimeA251Snake<-computeThermalizationTime(thermalizationDataA251Snake)
# rm(thermalizationDataA251Snake)
# 
# thermalizationDataB51Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=51_A=2.txt")
# thermalizationTimeB51Snake<-computeThermalizationTime(thermalizationDataB51Snake)
# rm(thermalizationDataB51Snake)
# thermalizationDataB101Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=101_A=2.txt")
# thermalizationTimeB101Snake<-computeThermalizationTime(thermalizationDataB101Snake)
# rm(thermalizationDataB101Snake)
# thermalizationDataB151Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=151_A=2.txt")
# thermalizationTimeB151Snake<-computeThermalizationTime(thermalizationDataB151Snake)
# rm(thermalizationDataB151Snake)
# thermalizationDataB201Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=201_A=2.txt")
# thermalizationTimeB201Snake<-computeThermalizationTime(thermalizationDataB201Snake)
# rm(thermalizationDataB201Snake)
# thermalizationDataB251Snake<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=251_A=2.txt")
# thermalizationTimeB251Snake<-computeThermalizationTime(thermalizationDataB251Snake)
# rm(thermalizationDataB251Snake)
# #####
# 
# #Pivot thermalization data
# ##########################
# thermalizationDataA51Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=51_A=3.txt")
# thermalizationTimeA51Pivot<-computeThermalizationTime(thermalizationDataA51Pivot)
# rm(thermalizationDataA51Pivot)
# thermalizationDataA101Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=101_A=3.txt")
# thermalizationTimeA101Pivot<-computeThermalizationTime(thermalizationDataA101Pivot)
# rm(thermalizationDataA101Pivot)
# thermalizationDataA151Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=151_A=3.txt")
# thermalizationTimeA151Pivot<-computeThermalizationTime(thermalizationDataA151Pivot)
# rm(thermalizationDataA151Pivot)
# thermalizationDataA201Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=201_A=3.txt")
# thermalizationTimeA201Pivot<-computeThermalizationTime(thermalizationDataA201Pivot)
# rm(thermalizationDataA201Pivot)
# thermalizationDataA251Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partA_S=251_A=3.txt")
# thermalizationTimeA251Pivot<-computeThermalizationTime(thermalizationDataA251Pivot)
# rm(thermalizationDataA251Pivot)
# 
# thermalizationDataB51Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=51_A=3.txt")
# thermalizationTimeB51Pivot<-computeThermalizationTime(thermalizationDataB51Pivot)
# rm(thermalizationDataB51Pivot)
# thermalizationDataB101Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=101_A=3.txt")
# thermalizationTimeB101Pivot<-computeThermalizationTime(thermalizationDataB101Pivot)
# rm(thermalizationDataB101Pivot)
# thermalizationDataB151Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=151_A=3.txt")
# thermalizationTimeB151Pivot<-computeThermalizationTime(thermalizationDataB151Pivot)
# rm(thermalizationDataB151Pivot)
# thermalizationDataB201Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=201_A=3.txt")
# thermalizationTimeB201Pivot<-computeThermalizationTime(thermalizationDataB201Pivot)
# rm(thermalizationDataB201Pivot)
# thermalizationDataB251Pivot<-read.csv("~/Documents/Huiswerk/Computational Physics/Random Walks/Data/thermalization_data_partB_S=251_A=3.txt")
# thermalizationTimeB251Pivot<-computeThermalizationTime(thermalizationDataB251Pivot)
# rm(thermalizationDataB251Pivot)
# #####
# 
# #Cat thermalization times
# #########################
# cat(thermalizationTimeA51Snake, file="Thermalization_Time_partA_S=51_A=Snake.txt", sep="\n")
# cat(thermalizationTimeA101Snake, file="Thermalization_Time_partA_S=101_A=Snake.txt", sep="\n")
# cat(thermalizationTimeA151Snake, file="Thermalization_Time_partA_S=151_A=Snake.txt", sep="\n")
# cat(thermalizationTimeA201Snake, file="Thermalization_Time_partA_S=201_A=Snake.txt", sep="\n")
# cat(thermalizationTimeA251Snake, file="Thermalization_Time_partA_S=251_A=Snake.txt", sep="\n")
# 
# cat(thermalizationTimeB51Snake, file="Thermalization_Time_partB_S=51_A=Snake.txt", sep="\n")
# cat(thermalizationTimeB101Snake, file="Thermalization_Time_partB_S=101_A=Snake.txt", sep="\n")
# cat(thermalizationTimeB151Snake, file="Thermalization_Time_partB_S=151_A=Snake.txt", sep="\n")
# cat(thermalizationTimeB201Snake, file="Thermalization_Time_partB_S=201_A=Snake.txt", sep="\n")
# cat(thermalizationTimeB251Snake, file="Thermalization_Time_partB_S=251_A=Snake.txt", sep="\n")
# 
# cat(thermalizationTimeA51Pivot, file="Thermalization_Time_partA_S=51_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeA101Pivot, file="Thermalization_Time_partA_S=101_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeA151Pivot, file="Thermalization_Time_partA_S=151_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeA201Pivot, file="Thermalization_Time_partA_S=201_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeA251Pivot, file="Thermalization_Time_partA_S=251_A=Pivot.txt", sep="\n")
# 
# cat(thermalizationTimeB51Pivot, file="Thermalization_Time_partB_S=51_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeB101Pivot, file="Thermalization_Time_partB_S=101_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeB151Pivot, file="Thermalization_Time_partB_S=151_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeB201Pivot, file="Thermalization_Time_partB_S=201_A=Pivot.txt", sep="\n")
# cat(thermalizationTimeB251Pivot, file="Thermalization_Time_partB_S=251_A=Pivot.txt", sep="\n")
# #####

#Read thermalization times
##########################
#thermalizationTimes<-readThermalizationTimes()
#####

#Local move algorithm
#################################################

#Read local data
usedSARWLengthsLocal<-10*c(1:7)+1
dataLocal<-readData(algorithmChoice=1)

#Compute autocorrelation function
autocorrelationEndToEndLocal<-computeEndToEndAutocorrelation(dataLocal[
  which(dataLocal$Size==max(usedSARWLengthsLocal)),])

#Guess the linear part of the correlation function
linearPartCorrelationLocal<-4500

#Plot the autocorrelation function with linear fit
autocorrelationEndToEndLocal.df<-data.frame(y=autocorrelationEndToEndLocal[1:min(which(autocorrelationEndToEndLocal<=0)-1)],
                                            x=seq(min(which(autocorrelationEndToEndLocal<=0)-1)))
linearFitCoefficientsLocal<-as.vector(coef(lm(log(autocorrelationEndToEndLocal.df[,1])~autocorrelationEndToEndLocal.df[,2], 
                   subset=c(1:linearPartCorrelationLocal))))

plotLocalEndToEndLengthCorrelation<-ggplot(autocorrelationEndToEndLocal.df, aes(x=x, y=y))+
  geom_point(aes(colour="correlatation"))+
  geom_line(aes(y=exp(linearFitCoefficientsLocal[2]*x+linearFitCoefficientsLocal[1]), colour="linear fit"))+
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(paste(Delta,italic(t))), y=expression(italic(C)[r[e]^2]))+
  theme_linedraw()
print(plotLocalEndToEndLengthCorrelation)
jpeg("./localEndToEndLengthCorrelationN=70.jpg")
print(plotLocalEndToEndLengthCorrelation)
dev.off()

#Compute correlation time
correlationTimeLocal<--1/linearFitCoefficientsLocal[2]

#Compute average end-to-end lengths
averageEndToEndLengthLocal<-computeAverageEndToEndLength(dataLocal, usedSARWLengthsLocal)
 
# #Plot average end-to-end lengths
#####
# plotLocalAverageEndToEndLength<-ggplot()+
#   geom_point(aes(y=averageEndToEndLengthLocal, x=usedSARWLengthsLocal))+
#   labs(x=expression(italic(N)), y=expression(paste("<",italic(r)[italic(e)]^2,"(",italic(N),")>")))+
#   scale_colour_manual(values=myPalette, guide=FALSE)+
#   theme_linedraw()
# print(plotLocalAverageEndToEndLength)
# jpeg("./localEndToEndLengthAverage.jpg")
# print(plotLocalAverageEndToEndLength)
# dev.off()
#####

#Compute growth exponent, nu
coefficientsEndToEndLocal<-computeExponent(averageEndToEndLengthLocal, usedSARWLengthsLocal)
nuLocal<-coefficientsEndToEndLocal[2]/2


#Plot average end-to-end lengths on log-log scale with linear fit
plotLocalEndToEndLengthFit<-ggplot()+
  geom_point(aes(x=usedSARWLengthsLocal-1, y=averageEndToEndLengthLocal, colour="averages"))+
  geom_line(aes(x=usedSARWLengthsLocal-1, y=exp(coefficientsEndToEndLocal[2]*log(usedSARWLengthsLocal)
                                                +coefficientsEndToEndLocal[1]), 
            colour="linear fit"))+
  scale_x_log10(breaks=seq(10, 70, 10)) + scale_y_log10(breaks=c(30,100,200,300,400,500,600))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(italic(N)), y=expression(paste("<",italic(r)[italic(e)]^2,">")))+
  theme_linedraw()
print(plotLocalEndToEndLengthFit)
jpeg("./localEndToEndLengthFit.jpg")
print(plotLocalEndToEndLengthFit)
dev.off()


rm(dataLocal)
#Slithering snake algorithm
#################################################

#Read snake data
usedSARWLengthsSnake<-10*c(1:15)+1
dataSnake<-readData(algorithmChoice=2)

#Compute autocorrelation function
autocorrelationEndToEndSnake<-computeEndToEndAutocorrelation(dataSnake[
  which(dataSnake$Size==max(usedSARWLengthsSnake)),])

#Guess the linear part of the correlation function
linearPartCorrelationSnake<-19

#Plot the autocorrelation function with linear fit
autocorrelationEndToEndSnake.df<-data.frame(y=autocorrelationEndToEndSnake[1:min(which(autocorrelationEndToEndSnake<=0)-1)],
                                            x=seq(min(which(autocorrelationEndToEndSnake<=0)-1)))
linearFitCoefficientsSnake<-as.vector(coef(lm(log(autocorrelationEndToEndSnake.df[,1])~autocorrelationEndToEndSnake.df[,2], 
                                              subset=c(1:linearPartCorrelationSnake))))

plotSnakeEndToEndLengthCorrelation<-ggplot(autocorrelationEndToEndSnake.df, aes(x=x, y=y))+
  geom_point(aes(colour="correlatation"))+
  geom_line(aes(y=exp(linearFitCoefficientsSnake[2]*x+linearFitCoefficientsSnake[1]), colour="linear fit"))+
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(paste(Delta,italic(t))), y=expression(italic(C)[r[e]^2]))+
  theme_linedraw()
print(plotSnakeEndToEndLengthCorrelation)
jpeg("./snakeEndToEndLengthCorrelationN=150.jpg")
print(plotSnakeEndToEndLengthCorrelation)
dev.off()

#Compute correlation time
correlationTimeSnake<--1/linearFitCoefficientsSnake[2]

#Compute average end-to-end lengths
averageEndToEndLengthSnake<-computeAverageEndToEndLength(dataSnake, usedSARWLengthsSnake)

#Compute growth exponent, nu
coefficientsEndToEndSnake<-computeExponent(averageEndToEndLengthSnake, usedSARWLengthsSnake)
nuSnake<-coefficientsEndToEndSnake[2]/2

#Plot average end-to-end lengths on log-log scale with linear fit
plotSnakeEndToEndLengthFit<-ggplot()+
  geom_point(aes(x=usedSARWLengthsSnake-1, y=averageEndToEndLengthSnake, colour="averages"))+
  geom_line(aes(x=usedSARWLengthsSnake-1, y=exp(coefficientsEndToEndSnake[2]*log(usedSARWLengthsSnake)
                                                +coefficientsEndToEndSnake[1]), 
                colour="linear fit"))+
  scale_x_log10(breaks=seq(10, 150, 20)) + scale_y_log10(breaks=c(30,seq(100,1500,200)))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(italic(N)), y=expression(paste("<",italic(r)[italic(e)]^2,">")))+
  theme_linedraw()
print(plotSnakeEndToEndLengthFit)
jpeg("./snakeEndToEndLengthFit.jpg")
print(plotSnakeEndToEndLengthFit)
dev.off()

rm(dataSnake)


#Pivot algorithm
#################################################

#Read pivot data
usedSARWLengthsPivot<-union(union(10*c(1:15)+1, 50*c(1:20)+1), 10*c(96:99)+1)
dataPivot<-readData(algorithmChoice=3)

#Compute autocorrelation function
autocorrelationEndToEndPivot<-computeEndToEndAutocorrelation(dataPivot[
  which(dataPivot$Size==max(usedSARWLengthsPivot)),])

#Guess the linear part of the correlation function
linearPartCorrelationPivot<-2

#Plot the autocorrelation function with linear fit
autocorrelationEndToEndPivot.df<-data.frame(y=abs(autocorrelationEndToEndPivot[1:10]), x=seq(10))
linearFitCoefficientsPivot<-as.vector(coef(lm(log(autocorrelationEndToEndPivot.df[,1])~autocorrelationEndToEndPivot.df[,2], 
                                              subset=c(1:linearPartCorrelationPivot))))

plotPivotEndToEndLengthCorrelation<-ggplot(autocorrelationEndToEndPivot.df, aes(x=x, y=y))+
  geom_point(aes(colour="correlatation"))+
  geom_line(aes(y=exp(linearFitCoefficientsPivot[2]*x+linearFitCoefficientsPivot[1]), colour="linear fit"))+
  scale_y_log10(breaks=c(0.0001, 0.01, 0.1, 1), limit=c(0.0004,1))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(paste(Delta,italic(t))), y=expression(italic(C)[r[e]^2]))+
  theme_linedraw()
print(plotPivotEndToEndLengthCorrelation)
jpeg("./pivotEndToEndLengthCorrelationN=1000.jpg")
print(plotPivotEndToEndLengthCorrelation)
dev.off()

#Compute correlation time
correlationTimePivot<--1/linearFitCoefficientsPivot[2]

#Compute average end-to-end lengths
averageEndToEndLengthPivot<-computeAverageEndToEndLength(dataPivot, usedSARWLengthsPivot)

#Compute growth exponent, nu
coefficientsEndToEndPivot<-computeExponent(averageEndToEndLengthPivot, usedSARWLengthsPivot)
nuPivot<-coefficientsEndToEndPivot[2]/2

#Plot average end-to-end lengths on log-log scale with linear fit
plotPivotEndToEndLengthFit<-ggplot()+
  geom_point(aes(x=usedSARWLengthsPivot-1, y=averageEndToEndLengthPivot, colour="averages"))+
  geom_line(aes(x=usedSARWLengthsPivot-1, y=exp(coefficientsEndToEndPivot[2]*log(usedSARWLengthsPivot)
                                                +coefficientsEndToEndPivot[1]), 
                colour="linear fit"))+
  scale_x_log10(breaks=c(seq(10, 70, 20),100,seq(200,1000,200))) + scale_y_log10()+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(italic(N)), y=expression(paste("<",italic(r)[italic(e)]^2,">")))+
  theme_linedraw()
print(plotPivotEndToEndLengthFit)
jpeg("./pivotEndToEndLengthFit.jpg")
print(plotPivotEndToEndLengthFit)
dev.off()

#Compute average garation radii
averageGyrationRadiusPivot<-computeAverageGyrationRadius(dataPivot, usedSARWLengthsPivot)
ratioGyrationToSETELPivot<-averageGyrationRadiusPivot/averageEndToEndLengthPivot
KratioPivot<-ratioGyrationToSETELPivot[length(ratioGyrationToSETELPivot)]

#Plot the ratio between gyration radius and squared end-to-end length
plotPivotGyration<-ggplot()+
  geom_point(aes(x=usedSARWLengthsPivot-1, y=ratioGyrationToSETELPivot, colour="gyration radius"))+
  scale_y_continuous(limits=c(0,0.25))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(italic(N)), 
       y=expression(frac(paste("<",italic(r)[italic(g)]^2,">"),paste("<",italic(r)[italic(e)]^2,">"))))+
  theme_linedraw()
print(plotPivotGyration)
jpeg("./pivotRatioGyration.jpg")
print(plotPivotGyration)
dev.off()

rm(dataPivot)


#Pivot algorithm in three dimensions
#################################################

#Read pivot3D data
usedSARWLengths3D<-10*c(1:12)+1
data3D<-readData3D()

#Compute autocorrelation function
autocorrelationEndToEnd3D<-computeEndToEndAutocorrelation(data3D[which(data3D$Size==max(usedSARWLengths3D)),])

#Guess the linear part of the correlation function
linearPartCorrelation3D<-2

#Plot the autocorrelation function with linear fit
autocorrelationEndToEnd3D.df<-data.frame(y=abs(autocorrelationEndToEnd3D[1:10]), x=seq(10))
linearFitCoefficients3D<-as.vector(coef(lm(log(autocorrelationEndToEnd3D.df[,1])~autocorrelationEndToEnd3D.df[,2], 
                                              subset=c(1:linearPartCorrelation3D))))

plot3DEndToEndLengthCorrelation<-ggplot(autocorrelationEndToEnd3D.df, aes(x=x, y=y))+
  geom_point(aes(colour="correlatation"))+
  geom_line(aes(y=exp(linearFitCoefficients3D[2]*x+linearFitCoefficients3D[1]), colour="linear fit"))+
  scale_y_log10(breaks=c(0.0001, 0.01, 0.1, 1), limit=c(0.0004,1))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(paste(Delta,italic(t))), y=expression(italic(C)[r[e]^2]))+
  theme_linedraw()
print(plot3DEndToEndLengthCorrelation)
jpeg("./3DEndToEndLengthCorrelationN=1000.jpg")
print(plot3DEndToEndLengthCorrelation)
dev.off()

#Compute correlation time
correlationTime3D<--1/linearFitCoefficients3D[2]

#Compute average end-to-end lengths
averageEndToEndLength3D<-computeAverageEndToEndLength(data3D, usedSARWLengths3D)

#Compute growth exponent, nu
coefficientsEndToEnd3D<-computeExponent(averageEndToEndLength3D, usedSARWLengths3D)
nu3D<-coefficientsEndToEnd3D[2]/2

#Plot average end-to-end lengths on log-log scale with linear fit
plot3DEndToEndLengthFit<-ggplot()+
  geom_point(aes(x=usedSARWLengths3D-1, y=averageEndToEndLength3D, colour="averages"))+
  geom_line(aes(x=usedSARWLengths3D-1, y=exp(coefficientsEndToEnd3D[2]*log(usedSARWLengths3D)
                                                +coefficientsEndToEnd3D[1]), 
                colour="linear fit"))+
  scale_x_log10(breaks=c(seq(10, 70, 20),100,seq(200,1000,200))) + scale_y_log10(breaks=c(50,seq(100,500,100)))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(italic(N)), y=expression(paste("<",italic(r)[italic(e)]^2,">")))+
  theme_linedraw()
print(plot3DEndToEndLengthFit)
jpeg("./3DEndToEndLengthFit.jpg")
print(plot3DEndToEndLengthFit)
dev.off()

#Compute average garation radii
averageGyrationRadius3D<-computeAverageGyrationRadius(data3D, usedSARWLengths3D)
ratioGyrationToSETEL3D<-averageGyrationRadius3D/averageEndToEndLength3D
Kratio3D<-ratioGyrationToSETEL3D[length(ratioGyrationToSETEL3D)]


#Plot the ratio between gyration radius and squared end-to-end length
plot3DGyration<-ggplot()+
  geom_point(aes(x=usedSARWLengths3D-1, y=ratioGyrationToSETEL3D, colour="gyration radius"))+
  scale_y_continuous(limits=c(0,0.25))+
  scale_colour_manual(values=myPalette, guide=FALSE)+
  labs(x=expression(italic(N)), 
       y=expression(frac(paste("<",italic(r)[italic(g)]^2,">"),paste("<",italic(r)[italic(e)]^2,">"))))+
  theme_linedraw()
print(plot3DGyration)
jpeg("./3DRatioGyration.jpg")
print(plot3DGyration)
dev.off()

rm(data3D)



#Stop the clock!
time<-proc.time()-time
time



