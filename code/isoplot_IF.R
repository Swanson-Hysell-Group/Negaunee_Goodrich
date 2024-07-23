
setwd("Github/Vulcan_Iron_Formation/data/geochron/iolite")

library(dplyr)
library(IsoplotR)
library(scales)
library(stringr)
library(ggplot2)

##### Function for splitting up samples in output table from Iolite and getting into isoplotr format####

splt_data<-function(df,sample,prop){
  if(prop==TRUE){
    samp<-df%>%
      filter(str_detect(X,regex(sample)))%>%
      select(c("Final.U238.Pb206_mean","Final.U238.Pb206_prop2SE","Final.Pb207.Pb206_mean",
               "Final.Pb207.Pb206_prop2SE","rho.207Pb.206Pb.v.238U.206Pb")) %>%
      na.omit(samp)
    read.data(samp,method='U-Pb',format=2,ierr=2)
  }
  else{
    samp<-df%>%
      filter(str_detect(X,regex(sample)))%>%
      select(c("Final.U238.Pb206_mean","Final.U238.Pb206_2SE",
               "Final.Pb207.Pb206_mean","Final.Pb207.Pb206_2SE","rho.207Pb.206Pb.v.238U.206Pb")) %>%
      na.omit(samp)
    read.data(samp,method='U-Pb',format=2,ierr=2)
    
  }


}


ref_uncert<-function(df){
  df%>%
    mutate("full_six_eight"=((HFO_six_eight_2SE**2)+(Final.Pb206.U238_prop2SE**2))**0.5)%>%
  
   na.omit(data) 
}


# HFO Abs reference values which are taken as the reciprocal for the 238/06 
HFO_three_eight_mean=1/0.12216
HFO_three_eight_se=1/0.0021

HFO_seven_six_mean=0.9167
HFO_seven_six_se=0.0061



########### Read in iolite file #############

april_day_01<-read.csv('0424_day1.csv')
april_day_02<-read.csv('0424_day2.csv')
april_day_03<-read.csv('0424_day3.csv')
nov<-read.csv("nov.csv")
jk25_liam<-read.csv("jk25_liam.csv")

shared_cols<-intersect(colnames(april_day_01),colnames(april_day_02))

jk_cols<-intersect(colnames(april_day_02),colnames(nov))
all_data<-rbind(april_day_01[,shared_cols],april_day_02[,shared_cols],april_day_03[,shared_cols])

jk22_01<-splt_data(april_day_01,"JK22",TRUE)
jk22_02<-splt_data(april_day_02,"JK22",TRUE)
jk22_both<-splt_data(all_data,"JK22",TRUE)

jk23_01<-splt_data(april_day_01,"JK23",TRUE)
jk23_02<-splt_data(april_day_02,"JK23",TRUE)
jk23_both<-splt_data(all_data,"JK23",TRUE)


jk25_01<-splt_data(nov,"JK_25_G",TRUE)
jk25_02<-splt_data(april_day_02, "JK25_2",TRUE)
jk25_both<-splt_data(rbind(april_day_02[,jk_cols],nov[,jk_cols]),"JK25",TRUE)
#jk25_both<-splt_data(jk25_liam,"JK25",FALSE)


jk32<-splt_data(april_day_03,"JK32",TRUE)

jk34<-splt_data(april_day_03,"JK34",TRUE)

jk36<-splt_data(april_day_02,"JK36",TRUE)

jk10<-splt_data(april_day_03,"JK10",TRUE)
jk21<-splt_data(nov,"JK_21_G_Band",TRUE)
jk21_02<-splt_data(april_day_03,"JK21_band",TRUE)
vein<-splt_data(nov,"JK_21_G_Vein",TRUE)
vein_02<-splt_data(april_day_03,"JK21_vein",TRUE)
NF10<-splt_data(april_day_01,"NF10",TRUE)
NF9<-splt_data(april_day_01,"NF9",TRUE)
jk3<-splt_data(april_day_03,"JK3_mph",TRUE)
Fe_grab<-splt_data(nov,"Fe-Grab",TRUE)

#write.table(test,'/Users/anthonyfuentes/Desktop/test.csv',col.names = FALSE,row.names=FALSE,sep=",")

#write.table(test_v2,'/Users/anthonyfuentes/Desktop/test_v2.csv',col.names = FALSE,row.names=FALSE,sep=",")
#### Jaspilite Facies####

jk21_pts<-c(1,27,33,41,44,55,58,60)

#pdf("../../../figures/geochron/jk21_01.pdf")
concordia(jk21,type=2,tlim=(5300:1100),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,cex.lab=(2),cex.axis=(1.5))
#dev.off()
concordia(vein_02,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)

concordia(jk21_02,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)


concordia(jk3,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)
####NF####
NF10_pt<-c(7,10,11,17,6,2,45,51,24,56,23,
           26,46,25,47,52,50,54,49,20,5,33,21,18,30)

#pdf("../../../figures/geochron/NF10_01.pdf")
concordia(NF10,type=2,tlim=(5000:1200),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,hide=NF10_pt,common.Pb = 0,cex.lab=(2),cex.axis=(1.5))
#dev.off()

NF9_pts<-c(14,33,63,62,64,61,74,68,65,
           28,66,17,73,43,44,56,51,37,51,50)
#pdf("../../../figures/geochron/NF9_01.pdf")
concordia(NF9,type=2,(5000:1200),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=NF9_pts,cex.lab=(2),cex.axis=(1.5))
#dev.off()
###### JK22  ##########

# Try and color by spot size 


jk22_full<-concordia(jk22_01,type=2,ellipse.fill = alpha('#de2d26',0.3),ticks = 10,show.numbers = TRUE,
                     concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)

jk22_all<-concordia(jk22_both,type=2,ellipse.fill = alpha('#de2d26',0.3),ticks = 10,show.numbers = TRUE,
                     concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)
####JK23 plotting by day first####

# First plot the full data set then second day. Finally, concatenate and plot



#pdf("../../../figures/geochron/jk23_full.pdf")
concordia(jk23_01,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)
#dev.off()

#pdf("../../../figures/geochron/jk23.pdf");
pts<-c(4,51,6,22,29,38,61,53,60,62,10,68,47,25,
       44,7,49,46,50,36,2,63,23,13,57)
concordia(jk23_01,type=2,tlim=(5300:540),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,omit=pts,cex.lab=(2),cex.axis=(1.5))
#dev.off();


#pdf("../../../figures/geochron/jk23day2_full.pdf");
concordia(jk23_02,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers =TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)
#dev.off()


#pdf("../../../figures/geochron/jk23day2_curated.pdf");
pts_jk02<-c(59,60,32,54,46,45,55,2,8,21,3,58,18,13,
            15,34,56,35,20,37,23,48,36,30,29,47,12,
            51,26,28,22,52,50,53,9,31)
concordia(jk23_02,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=pts_jk02)
#dev.off()

pts_jk23<-c(12,123,113,13,38,29,82,22,102,83,127,4,98,70,63,
            76,126,103,97,89,53,62,10,46,100,12,60,115,117,116,104,2,96,90,
            85,92,94,75,51,25,49,128,36,86,105,6,106,8,112,
            47,114,122,124,71,91,7,35,69,87,61,50,120)

#pdf("../../../figures/geochron/jk23_final.pdf")
concordia(jk23_both,type=2,tlim=(5000:500),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=pts_jk23,cex.lab=(2),cex.axis=(1.5))
#dev.off()


####JK 25####


# Once I have Liam's iolite file I can export the second day data and merge values based on shared values 
#to pull out the excess uncertainty parameter 
pdf("../../../figures/geochron/jk25.pdf")
concordia(jk25_01,type=2,tlim=(5300:550),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers =FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,omit=c(27,1,46),cex.lab=(2),cex.axis=(1.5))
dev.off()
concordia(jk25_02,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers=TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=)


concordia(jk25_both,type=2,tlim=(5300:600),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers =FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=c(77),cex.lab=(2),cex.axis=(1.5))


#### Other JK samples####
pts_jk32<-c(1,3,4,9,10,18,25,26,27,33,21,43,44,42,45,47,48,39,41,46,49,50,52)

concordia(jk32,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=c(45,44,42,43))
pdf("../../../figures/geochron/jk32.pdf")
concordia(jk32,type=2,tlim=(5300:600),,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=pts_jk32,cex.lab=(2),cex.axis=(1.5))
dev.off()


concordia(jk34,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=c(9,33,30,20,28,2,29,37))

concordia(jk36,type=2,ellipse.fill = alpha('#a6bddb',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)




#### Goodrich Conglomerate####

GQ6_01<-splt_data(april_day_01,"GQ6_2_matrix",TRUE)

#pts_gq1<-c(13,3,17)
pdf("../../../figures/geochron/GQ6_2_day1.pdf")
concordia(GQ6_01,type=2,tlim=(5300:1500),ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers = FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,cex.lab=(2),cex.axis=(1.5))

dev.off()



GQ6_02<-splt_data(april_day_02,"GQ6_2_matrix",TRUE)

pts_GQ_01<-c(55,54,51,50,23,26,53,48,49,27,52,2,3,4,
          22,24,21,25,5,6,1,20,38)
pts_GQ<-c(31,32,12,11,2,4,7,2,28,45,47,18,36,20,23,28,44,46,10,13,
          43,9,8,29,30,42,31,32,33,34,35,41,35,40,39,38,37,19,17,16,15,14)
concordia(GQ6_02,type=2,tlim=(5200:1850),ellipse.fill = alpha('#a6bddb',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=pts_GQ_01)

GQ_mrt<-splt_data(april_day_01,"GQ6_2_mrt",TRUE)
concordia(GQ_mrt,type=2,ellipse.fill = alpha('#a6bddb',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0)


GQ_both<-splt_data(all_data,"GQ6_2_matrix",TRUE)

gq_both_pts<-c(66,67,107,93,73,68,65,110,109,106,
               105,81,63,64,101,99,101,78,79,104,55,77,108,98,
               77)

concordia(GQ_both,type=2,ellipse.fill = alpha('#a6bddb',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=gq_both_pts)


GQ2<-splt_data(nov,"GQ2",TRUE)

GQ2_pts<-c(9,26,12,30,2,23,16,14,5,35,34,21,1,15,28,19,32,18,
           27,13,6,33,20,10,31)
concordia(GQ2,type=2,tlim=(5200:1650),ellipse.fill = alpha('#a6bddb',0.6),ticks = 10,show.numbers = TRUE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,hide=GQ2_pts)


#### Fe Grab####

fe_pts<-c(1,16,36,18,4,47,44,46,42,66,43,
          20,22,29,45,6,31,34,39,12,37,40,7,12,24)

#pdf("../../../figures/geochron/Fe_Grab.pdf")
concordia(Fe_grab,type=2,ellipse.fill = alpha('#de2d26',0.6),ticks = 10,show.numbers =FALSE,
          concordia.col = 'black',show.age = 2,oerr=2,common.Pb = 0,omit=fe_pts,,cex.lab=(2),cex.axis=(1.5))
#dev.off()
