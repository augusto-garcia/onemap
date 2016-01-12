### R code from vignette source 'Tutorial_Onemap_reduced_version.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Tutorial_Onemap_reduced_version.Rnw:93-94 (eval = FALSE)
###################################################
## options(width=40, show.signif.stars=FALSE)


###################################################
### code chunk number 2: Tutorial_Onemap_reduced_version.Rnw:98-99 (eval = FALSE)
###################################################
## 2+3


###################################################
### code chunk number 3: Tutorial_Onemap_reduced_version.Rnw:102-103 (eval = FALSE)
###################################################
## x<-2+3


###################################################
### code chunk number 4: Tutorial_Onemap_reduced_version.Rnw:106-107 (eval = FALSE)
###################################################
## x


###################################################
### code chunk number 5: Tutorial_Onemap_reduced_version.Rnw:110-111 (eval = FALSE)
###################################################
## x+4


###################################################
### code chunk number 6: Tutorial_Onemap_reduced_version.Rnw:117-118 (eval = FALSE)
###################################################
## log(6.7)


###################################################
### code chunk number 7: Tutorial_Onemap_reduced_version.Rnw:121-122 (eval = FALSE)
###################################################
## log(6.7,base=4)


###################################################
### code chunk number 8: Tutorial_Onemap_reduced_version.Rnw:125-127 (eval = FALSE)
###################################################
## y<-c(6.7, 3.2, 5.4, 8.1, 4.9, 9.7, 2.5)
## log(y)


###################################################
### code chunk number 9: Tutorial_Onemap_reduced_version.Rnw:148-149 (eval = FALSE)
###################################################
## install.packages("onemap")


###################################################
### code chunk number 10: Tutorial_Onemap_reduced_version.Rnw:154-155 (eval = FALSE)
###################################################
## library("onemap")


###################################################
### code chunk number 11: Tutorial_Onemap_reduced_version.Rnw:160-161 (eval = FALSE)
###################################################
## rf<-c(0.01, 0.12, 0.05, 0.11, 0.21, 0.07)


###################################################
### code chunk number 12: Tutorial_Onemap_reduced_version.Rnw:166-167 (eval = FALSE)
###################################################
## kosambi(rf)


###################################################
### code chunk number 13: Tutorial_Onemap_reduced_version.Rnw:202-203 (eval = FALSE)
###################################################
## setwd("C:/Users/mmollina/Documents")


###################################################
### code chunk number 14: Tutorial_Onemap_reduced_version.Rnw:210-211 (eval = FALSE)
###################################################
## (dat<-read.table(file="test.txt", header=TRUE))


###################################################
### code chunk number 15: Tutorial_Onemap_reduced_version.Rnw:214-218 (eval = FALSE)
###################################################
## dat<-data.frame(
## x=c(2.13,4.48,10.95,10.03,12.72,24.63,22.57,29.78,19.54,7.86,11.75,23.71),
## y=c(4.50,1.98,9.29,16.25,27.38,22.60,36.87,31.73,10.42,14.68,8.68,37.39))
## dat


###################################################
### code chunk number 16: Tutorial_Onemap_reduced_version.Rnw:223-225 (eval = FALSE)
###################################################
## dat$x
## dat$y 


###################################################
### code chunk number 17: Tutorial_Onemap_reduced_version.Rnw:230-233 (eval = FALSE)
###################################################
## summary(dat)
## summary(dat$x)
## summary(dat$y)


###################################################
### code chunk number 18: Tutorial_Onemap_reduced_version.Rnw:238-239 (eval = FALSE)
###################################################
## write.table(x=summary(dat), file="test_sum.txt", quote=FALSE)


###################################################
### code chunk number 19: Tutorial_Onemap_reduced_version.Rnw:248-249 (eval = FALSE)
###################################################
## class(dat)


###################################################
### code chunk number 20: Tutorial_Onemap_reduced_version.Rnw:254-256 (eval = FALSE)
###################################################
## ft.mod<-lm(dat$y~dat$x)
## ft.mod


###################################################
### code chunk number 21: Tutorial_Onemap_reduced_version.Rnw:261-262 (eval = FALSE)
###################################################
## class(ft.mod)


###################################################
### code chunk number 22: Tutorial_Onemap_reduced_version.Rnw:267-268 (eval = FALSE)
###################################################
## summary(ft.mod)


###################################################
### code chunk number 23: Tutorial_Onemap_reduced_version.Rnw:279-280 (eval = FALSE)
###################################################
## save.image("myworkspace.RData")


###################################################
### code chunk number 24: Tutorial_Onemap_reduced_version.Rnw:285-286 (eval = FALSE)
###################################################
## load("myworkspace.RData")


###################################################
### code chunk number 25: Tutorial_Onemap_reduced_version.Rnw:295-296 (eval = FALSE)
###################################################
## install.packages("onemap")


###################################################
### code chunk number 26: Tutorial_Onemap_reduced_version.Rnw:352-353 (eval = FALSE)
###################################################
## library(onemap)


###################################################
### code chunk number 27: Tutorial_Onemap_reduced_version.Rnw:358-360 (eval = FALSE)
###################################################
## library()
## data()


###################################################
### code chunk number 28: Tutorial_Onemap_reduced_version.Rnw:371-372 (eval = FALSE)
###################################################
## library(onemap)


###################################################
### code chunk number 29: Tutorial_Onemap_reduced_version.Rnw:376-377 (eval = FALSE)
###################################################
## save.image("C:/.../yourfile.RData")


###################################################
### code chunk number 30: Tutorial_Onemap_reduced_version.Rnw:445-446 (eval = FALSE)
###################################################
## system.file(package="onemap")


###################################################
### code chunk number 31: Tutorial_Onemap_reduced_version.Rnw:452-453 (eval = FALSE)
###################################################
## example.out<- read.outcross("C:/workingdirectory","example.out.txt")


###################################################
### code chunk number 32: Tutorial_Onemap_reduced_version.Rnw:455-456 (eval = FALSE)
###################################################
## example.out<- read.outcross(system.file("example",package="onemap"),"example.out.txt")


###################################################
### code chunk number 33: Tutorial_Onemap_reduced_version.Rnw:461-462 (eval = FALSE)
###################################################
## example.out<- read.outcross(file="example.out.txt")


###################################################
### code chunk number 34: Tutorial_Onemap_reduced_version.Rnw:466-467 (eval = FALSE)
###################################################
## example.out<- read.outcross(system.file("example",package="onemap"),"example.out.txt")


###################################################
### code chunk number 35: Tutorial_Onemap_reduced_version.Rnw:471-472 (eval = FALSE)
###################################################
## data(example.out)


###################################################
### code chunk number 36: Tutorial_Onemap_reduced_version.Rnw:476-477 (eval = FALSE)
###################################################
## example.out


###################################################
### code chunk number 37: Tutorial_Onemap_reduced_version.Rnw:481-482 (eval = FALSE)
###################################################
## example.out


###################################################
### code chunk number 38: Tutorial_Onemap_reduced_version.Rnw:490-491 (eval = FALSE)
###################################################
## twopts <- rf.2pts(example.out)


###################################################
### code chunk number 39: Tutorial_Onemap_reduced_version.Rnw:496-497 (eval = FALSE)
###################################################
## twopts <- rf.2pts(example.out, LOD=3, max.rf=0.4)


###################################################
### code chunk number 40: Tutorial_Onemap_reduced_version.Rnw:503-504 (eval = FALSE)
###################################################
## twopts


###################################################
### code chunk number 41: Tutorial_Onemap_reduced_version.Rnw:508-509 (eval = FALSE)
###################################################
## twopts


###################################################
### code chunk number 42: Tutorial_Onemap_reduced_version.Rnw:513-514 (eval = FALSE)
###################################################
## print(twopts, "M1", "M3")


###################################################
### code chunk number 43: Tutorial_Onemap_reduced_version.Rnw:524-525 (eval = FALSE)
###################################################
## mark.all <- make.seq(twopts, "all")


###################################################
### code chunk number 44: Tutorial_Onemap_reduced_version.Rnw:529-530 (eval = FALSE)
###################################################
## marker.type(mark.all)


###################################################
### code chunk number 45: Tutorial_Onemap_reduced_version.Rnw:534-535 (eval = FALSE)
###################################################
## LGs <- group(mark.all)


###################################################
### code chunk number 46: Tutorial_Onemap_reduced_version.Rnw:540-541 (eval = FALSE)
###################################################
## LGs


###################################################
### code chunk number 47: Tutorial_Onemap_reduced_version.Rnw:545-546 (eval = FALSE)
###################################################
## LGs


###################################################
### code chunk number 48: Tutorial_Onemap_reduced_version.Rnw:550-551 (eval = FALSE)
###################################################
## print(LGs, detailed=FALSE)


###################################################
### code chunk number 49: Tutorial_Onemap_reduced_version.Rnw:555-557 (eval = FALSE)
###################################################
## LGs <- group(mark.all, LOD=6)
## LGs


###################################################
### code chunk number 50: Tutorial_Onemap_reduced_version.Rnw:561-563 (eval = FALSE)
###################################################
## LGs <- group(mark.all, LOD=3, max.rf=0.4)
## LGs


###################################################
### code chunk number 51: Tutorial_Onemap_reduced_version.Rnw:571-572 (eval = FALSE)
###################################################
## set.map.fun(type="haldane")


###################################################
### code chunk number 52: Tutorial_Onemap_reduced_version.Rnw:576-577 (eval = FALSE)
###################################################
## set.map.fun(type="kosambi")


###################################################
### code chunk number 53: Tutorial_Onemap_reduced_version.Rnw:581-582 (eval = FALSE)
###################################################
## LG3 <- make.seq(LGs, 3)


###################################################
### code chunk number 54: Tutorial_Onemap_reduced_version.Rnw:587-588 (eval = FALSE)
###################################################
## LG3


###################################################
### code chunk number 55: Tutorial_Onemap_reduced_version.Rnw:591-592 (eval = FALSE)
###################################################
## LG3


###################################################
### code chunk number 56: Tutorial_Onemap_reduced_version.Rnw:597-601 (eval = FALSE)
###################################################
## LG3.ser <- seriation(LG3)
## LG3.rcd <- rcd(LG3)
## LG3.rec <- record(LG3)
## LG3.ug  <- ug(LG3)


###################################################
### code chunk number 57: Tutorial_Onemap_reduced_version.Rnw:607-608 (eval = FALSE)
###################################################
## LG3.comp <- compare(LG3)


###################################################
### code chunk number 58: Tutorial_Onemap_reduced_version.Rnw:616-617 (eval = FALSE)
###################################################
## LG3.comp


###################################################
### code chunk number 59: Tutorial_Onemap_reduced_version.Rnw:627-628 (eval = FALSE)
###################################################
## LG3.final <- make.seq(LG3.comp,1,1)


###################################################
### code chunk number 60: Tutorial_Onemap_reduced_version.Rnw:634-635 (eval = FALSE)
###################################################
## LG3.final <- make.seq(LG3.comp)


###################################################
### code chunk number 61: Tutorial_Onemap_reduced_version.Rnw:640-641 (eval = FALSE)
###################################################
## LG3.final


###################################################
### code chunk number 62: Tutorial_Onemap_reduced_version.Rnw:653-655 (eval = FALSE)
###################################################
## LG2 <- make.seq(LGs, 2)
## LG2


###################################################
### code chunk number 63: Tutorial_Onemap_reduced_version.Rnw:660-662 (eval = FALSE)
###################################################
## LG2.rcd <- rcd(LG2)
## LG2.rcd


###################################################
### code chunk number 64: Tutorial_Onemap_reduced_version.Rnw:666-667 (eval = FALSE)
###################################################
## marker.type(LG2)


###################################################
### code chunk number 65: Tutorial_Onemap_reduced_version.Rnw:671-673 (eval = FALSE)
###################################################
## LG2.init <- make.seq(twopts,c(4,23,19,20,24))
## LG2.comp <- compare(LG2.init)


###################################################
### code chunk number 66: Tutorial_Onemap_reduced_version.Rnw:675-676 (eval = FALSE)
###################################################
## LG2.comp


###################################################
### code chunk number 67: Tutorial_Onemap_reduced_version.Rnw:681-682 (eval = FALSE)
###################################################
## LG2.frame <- make.seq(LG2.comp)


###################################################
### code chunk number 68: Tutorial_Onemap_reduced_version.Rnw:686-688 (eval = FALSE)
###################################################
## LG2.extend <- try_seq(LG2.frame,9)
## LG2.extend


###################################################
### code chunk number 69: Tutorial_Onemap_reduced_version.Rnw:693-694 (eval = FALSE)
###################################################
## print(LG2.extend,5)


###################################################
### code chunk number 70: Tutorial_Onemap_reduced_version.Rnw:701-702 (eval = FALSE)
###################################################
## LG2.extend <- try_seq(LG2.frame, 9, draw.try=TRUE)


###################################################
### code chunk number 71: Tutorial_Onemap_reduced_version.Rnw:708-709 (eval = FALSE)
###################################################
## LG2.frame <- make.seq(LG2.extend,5,1)


###################################################
### code chunk number 72: Tutorial_Onemap_reduced_version.Rnw:717-725 (eval = FALSE)
###################################################
## LG2.extend <- try_seq(LG2.frame,29)
## LG2.frame <- make.seq(LG2.extend,7)
## LG2.extend <- try_seq(LG2.frame,27)
## LG2.frame <- make.seq(LG2.extend,1)
## LG2.extend <- try_seq(LG2.frame, 16)
## LG2.frame <- make.seq(LG2.extend,2)
## LG2.extend <- try_seq(LG2.frame,21)
## LG2.final <- make.seq(LG2.extend,6)


###################################################
### code chunk number 73: Tutorial_Onemap_reduced_version.Rnw:730-731 (eval = FALSE)
###################################################
## LG2.ord <- order.seq(LG2, n.init=5, THRES=3, draw.try=TRUE, wait=1)


###################################################
### code chunk number 74: Tutorial_Onemap_reduced_version.Rnw:739-740 (eval = FALSE)
###################################################
## LG2.ord


###################################################
### code chunk number 75: Tutorial_Onemap_reduced_version.Rnw:745-746 (eval = FALSE)
###################################################
## LG2.safe <- make.seq(LG2.ord,"safe")


###################################################
### code chunk number 76: Tutorial_Onemap_reduced_version.Rnw:749-751 (eval = FALSE)
###################################################
## LG2.all <- make.seq(LG2.ord,"force")
## LG2.all


###################################################
### code chunk number 77: Tutorial_Onemap_reduced_version.Rnw:756-758 (eval = FALSE)
###################################################
## LG2.ord <- order.seq(LG2, n.init=5, THRES=3, touchdown=TRUE)
## LG2.ord


###################################################
### code chunk number 78: Tutorial_Onemap_reduced_version.Rnw:763-764 (eval = FALSE)
###################################################
## ripple_seq(LG2.all, ws=4, LOD=3)


###################################################
### code chunk number 79: Tutorial_Onemap_reduced_version.Rnw:775-776 (eval = FALSE)
###################################################
## LG2.all


###################################################
### code chunk number 80: Tutorial_Onemap_reduced_version.Rnw:784-785 (eval = FALSE)
###################################################
## LG1 <- make.seq(LGs, 1)


###################################################
### code chunk number 81: Tutorial_Onemap_reduced_version.Rnw:789-791 (eval = FALSE)
###################################################
## LG1.ord <- order.seq(LG1, n.init=6, touchdown=TRUE)
## LG1.ord


###################################################
### code chunk number 82: Tutorial_Onemap_reduced_version.Rnw:796-797 (eval = FALSE)
###################################################
## (LG1.final <- make.seq(LG1.ord,"force"))


###################################################
### code chunk number 83: Tutorial_Onemap_reduced_version.Rnw:801-802 (eval = FALSE)
###################################################
## ripple_seq(LG1.final)


###################################################
### code chunk number 84: Tutorial_Onemap_reduced_version.Rnw:808-809 (eval = FALSE)
###################################################
## LG1.final


###################################################
### code chunk number 85: Tutorial_Onemap_reduced_version.Rnw:814-818 (eval = FALSE)
###################################################
## LG1.ser <- seriation(LG1)
## LG1.rcd <- rcd(LG1)
## LG1.rec <- record(LG1)
## LG1.ug  <- ug(LG1)


###################################################
### code chunk number 86: Tutorial_Onemap_reduced_version.Rnw:830-832 (eval = FALSE)
###################################################
## any.seq <- make.seq(twopts,c(30,12,3,14,2))
## (any.seq.map <- map(any.seq))


###################################################
### code chunk number 87: Tutorial_Onemap_reduced_version.Rnw:837-839 (eval = FALSE)
###################################################
## any.seq <- make.seq(twopts,c(30,12,3,14,2),phase=c(4,1,4,3))
## (any.seq.map <- map(any.seq))


###################################################
### code chunk number 88: Tutorial_Onemap_reduced_version.Rnw:843-844 (eval = FALSE)
###################################################
## (any.seq <- add.marker(any.seq, 4:8))


###################################################
### code chunk number 89: Tutorial_Onemap_reduced_version.Rnw:848-849 (eval = FALSE)
###################################################
## (any.seq <- drop.marker(any.seq, c(3,4,5,12,30)))


###################################################
### code chunk number 90: Tutorial_Onemap_reduced_version.Rnw:863-864 (eval = FALSE)
###################################################
## (LGs <- group(mark.all, LOD=2.5))


###################################################
### code chunk number 91: Tutorial_Onemap_reduced_version.Rnw:871-873 (eval = FALSE)
###################################################
## LG.err<-make.seq(LGs, 2)
## LG.err.ord<-order.seq(LG.err)


###################################################
### code chunk number 92: Tutorial_Onemap_reduced_version.Rnw:878-879 (eval = FALSE)
###################################################
## (LG.err.map<-make.seq(LG.err.ord, "force"))


###################################################
### code chunk number 93: Tutorial_Onemap_reduced_version.Rnw:882-883 (eval = FALSE)
###################################################
## (LG.err.map<-map(make.seq(twopts,c(27,16,20,4,19,21,23,9,24,29,22,7,18, 8, 13), phase=c(1,1,2,4,1,3,1,2,3,3,2,1,4,4))))


###################################################
### code chunk number 94: Tutorial_Onemap_reduced_version.Rnw:888-889 (eval = FALSE)
###################################################
## rf.graph.table(LG.err.map)


###################################################
### code chunk number 95: Tutorial_Onemap_reduced_version.Rnw:892-893 (eval = FALSE)
###################################################
## rf.graph.table(LG.err.map, inter=FALSE)


###################################################
### code chunk number 96: Tutorial_Onemap_reduced_version.Rnw:896-897 (eval = FALSE)
###################################################
## dev.off()


###################################################
### code chunk number 97: Tutorial_Onemap_reduced_version.Rnw:916-918 (eval = FALSE)
###################################################
## maps<-list(LG1.final, LG2.final, LG3.final)
## draw.map(maps, names= TRUE, grid=TRUE, cex.mrk=0.7)


###################################################
### code chunk number 98: Tutorial_Onemap_reduced_version.Rnw:923-924 (eval = FALSE)
###################################################
## draw.map(LG1.final, names= TRUE, grid=TRUE, cex.mrk=0.7)


###################################################
### code chunk number 99: Tutorial_Onemap_reduced_version.Rnw:999-1000 (eval = FALSE)
###################################################
## library(onemap)


###################################################
### code chunk number 100: Tutorial_Onemap_reduced_version.Rnw:1004-1005 (eval = FALSE)
###################################################
## save.image("C:/.../yourfile.RData")


###################################################
### code chunk number 101: Tutorial_Onemap_reduced_version.Rnw:1017-1019 (eval = FALSE)
###################################################
## fake.f2.onemap <- read.mapmaker(dir="C:/workingdirectory", 
##                                 file="your_data_file.raw")


###################################################
### code chunk number 102: Tutorial_Onemap_reduced_version.Rnw:1026-1028 (eval = FALSE)
###################################################
## data(fake.f2.onemap)
## fake.f2.onemap


###################################################
### code chunk number 103: Tutorial_Onemap_reduced_version.Rnw:1039-1040 (eval = FALSE)
###################################################
## twopts.f2 <- rf.2pts(fake.f2.onemap)


###################################################
### code chunk number 104: Tutorial_Onemap_reduced_version.Rnw:1046-1047 (eval = FALSE)
###################################################
## print(twopts.f2, "M12", "M42")


###################################################
### code chunk number 105: Tutorial_Onemap_reduced_version.Rnw:1057-1058 (eval = FALSE)
###################################################
## mark.all.f2 <- make.seq(twopts.f2, "all")


###################################################
### code chunk number 106: Tutorial_Onemap_reduced_version.Rnw:1062-1063 (eval = FALSE)
###################################################
## mrk.subset<-make.seq(twopts.f2, c(1,3,7))


###################################################
### code chunk number 107: Tutorial_Onemap_reduced_version.Rnw:1067-1068 (eval = FALSE)
###################################################
## (LGs.f2 <- group(mark.all.f2, LOD=3, max.rf=0.5))


###################################################
### code chunk number 108: Tutorial_Onemap_reduced_version.Rnw:1079-1080 (eval = FALSE)
###################################################
## set.map.fun(type="haldane")


###################################################
### code chunk number 109: Tutorial_Onemap_reduced_version.Rnw:1084-1085 (eval = FALSE)
###################################################
## set.map.fun(type="kosambi")


###################################################
### code chunk number 110: Tutorial_Onemap_reduced_version.Rnw:1090-1091 (eval = FALSE)
###################################################
## LG2.f2 <- make.seq(LGs.f2, 2)


###################################################
### code chunk number 111: Tutorial_Onemap_reduced_version.Rnw:1097-1098 (eval = FALSE)
###################################################
## LG2.f2


###################################################
### code chunk number 112: Tutorial_Onemap_reduced_version.Rnw:1101-1102 (eval = FALSE)
###################################################
## LG2.f2


###################################################
### code chunk number 113: Tutorial_Onemap_reduced_version.Rnw:1107-1111 (eval = FALSE)
###################################################
## LG2.ser.f2 <- seriation(LG2.f2)
## LG2.rcd.f2 <- rcd(LG2.f2)
## LG2.rec.f2 <- record(LG2.f2)
## LG2.ug.f2 <- ug(LG2.f2)


###################################################
### code chunk number 114: Tutorial_Onemap_reduced_version.Rnw:1122-1127 (eval = FALSE)
###################################################
## LG2.f2.ord <- order.seq(input.seq=LG2.f2, n.init = 5, 
##                         subset.search = "twopt", 
##                         twopt.alg = "rcd", THRES = 3, 
##                         draw.try = TRUE, wait = 1)
##                         


###################################################
### code chunk number 115: Tutorial_Onemap_reduced_version.Rnw:1134-1135 (eval = FALSE)
###################################################
## LG2.f2.ord


###################################################
### code chunk number 116: Tutorial_Onemap_reduced_version.Rnw:1140-1141 (eval = FALSE)
###################################################
## LG2.f2.safe <- make.seq(LG2.f2.ord,"safe")


###################################################
### code chunk number 117: Tutorial_Onemap_reduced_version.Rnw:1144-1145 (eval = FALSE)
###################################################
## (LG2.f2.all <- make.seq(LG2.f2.ord,"force"))


###################################################
### code chunk number 118: Tutorial_Onemap_reduced_version.Rnw:1151-1156 (eval = FALSE)
###################################################
## LG2.f2.ord <- order.seq(input.seq=LG2.f2, n.init = 5, 
##                         subset.search = "twopt", 
##                         twopt.alg = "rcd", THRES = 3, 
##                         draw.try = TRUE, wait = 1,
##                         touchdown=TRUE)


###################################################
### code chunk number 119: Tutorial_Onemap_reduced_version.Rnw:1161-1162 (eval = FALSE)
###################################################
## (LG2.f2.final<-make.seq(LG2.f2.ord, "force"))


###################################################
### code chunk number 120: Tutorial_Onemap_reduced_version.Rnw:1166-1167 (eval = FALSE)
###################################################
## ripple_seq(LG2.f2.final, ws=5, LOD=3)


###################################################
### code chunk number 121: Tutorial_Onemap_reduced_version.Rnw:1177-1178 (eval = FALSE)
###################################################
## LG2.f2.final


###################################################
### code chunk number 122: Tutorial_Onemap_reduced_version.Rnw:1187-1188 (eval = FALSE)
###################################################
## LG1.f2 <- make.seq(LGs.f2, 1)


###################################################
### code chunk number 123: Tutorial_Onemap_reduced_version.Rnw:1192-1197 (eval = FALSE)
###################################################
## LG1.f2.ord <- order.seq(input.seq=LG1.f2, n.init = 5, 
##                         subset.search = "twopt", 
##                         twopt.alg = "rcd", THRES = 3, 
##                         draw.try = TRUE, wait = 1,
##                         touchdown=TRUE)


###################################################
### code chunk number 124: Tutorial_Onemap_reduced_version.Rnw:1202-1203 (eval = FALSE)
###################################################
## (LG1.f2.final <- make.seq(LG1.f2.ord,"force"))


###################################################
### code chunk number 125: Tutorial_Onemap_reduced_version.Rnw:1207-1208 (eval = FALSE)
###################################################
## ripple_seq(ws=5, LG1.f2.final)


###################################################
### code chunk number 126: Tutorial_Onemap_reduced_version.Rnw:1214-1215 (eval = FALSE)
###################################################
## LG1.f2.final


###################################################
### code chunk number 127: Tutorial_Onemap_reduced_version.Rnw:1225-1226 (eval = FALSE)
###################################################
## LG3.f2 <- make.seq(LGs.f2, 3)


###################################################
### code chunk number 128: Tutorial_Onemap_reduced_version.Rnw:1230-1235 (eval = FALSE)
###################################################
## LG3.f2.ord <- order.seq(input.seq=LG3.f2, n.init = 5, 
##                         subset.search = "twopt", 
##                         twopt.alg = "rcd", THRES = 3, 
##                         draw.try = TRUE, wait = 1,
##                         touchdown=TRUE)


###################################################
### code chunk number 129: Tutorial_Onemap_reduced_version.Rnw:1240-1241 (eval = FALSE)
###################################################
## (LG3.f2.final <- make.seq(LG3.f2.ord,"force"))


###################################################
### code chunk number 130: Tutorial_Onemap_reduced_version.Rnw:1245-1246 (eval = FALSE)
###################################################
## ripple_seq(ws=5, LG3.f2.final)


###################################################
### code chunk number 131: Tutorial_Onemap_reduced_version.Rnw:1252-1253 (eval = FALSE)
###################################################
## LG3.f2.final


###################################################
### code chunk number 132: Tutorial_Onemap_reduced_version.Rnw:1264-1266 (eval = FALSE)
###################################################
## LG3seq.f2 <- make.seq(twopts.f2,c(47,38,59,16,62,21,20,48,22))
## (LG3seq.f2.map <- map(LG3seq.f2))


###################################################
### code chunk number 133: Tutorial_Onemap_reduced_version.Rnw:1271-1272 (eval = FALSE)
###################################################
##  marker.type(LG3seq.f2.map)


###################################################
### code chunk number 134: Tutorial_Onemap_reduced_version.Rnw:1276-1277 (eval = FALSE)
###################################################
## (LG3seq.f2.map <- add.marker(LG3seq.f2.map, c(18,56,50)))


###################################################
### code chunk number 135: Tutorial_Onemap_reduced_version.Rnw:1281-1282 (eval = FALSE)
###################################################
## (LG3seq.f2.map <- drop.marker(LG3seq.f2.map, c(59,21)))


###################################################
### code chunk number 136: Tutorial_Onemap_reduced_version.Rnw:1293-1294 (eval = FALSE)
###################################################
## temp.seq<-drop.marker(LG3.f2.final, 38)


###################################################
### code chunk number 137: Tutorial_Onemap_reduced_version.Rnw:1297-1299 (eval = FALSE)
###################################################
## (temp.seq<-add.marker(temp.seq, 38))
## (LG3.f2.wrong<-map(temp.seq))


###################################################
### code chunk number 138: Tutorial_Onemap_reduced_version.Rnw:1306-1307 (eval = FALSE)
###################################################
## rf.graph.table(LG3.f2.wrong)


###################################################
### code chunk number 139: Tutorial_Onemap_reduced_version.Rnw:1310-1311 (eval = FALSE)
###################################################
## rf.graph.table(LG3.f2.wrong, inter=FALSE)


###################################################
### code chunk number 140: Tutorial_Onemap_reduced_version.Rnw:1320-1323 (eval = FALSE)
###################################################
## temp.seq <- drop.marker(LG3.f2.wrong,38)
## temp.map <- map(temp.seq)
## temp.try <- try_seq(temp.map, 38, draw.try=TRUE)


###################################################
### code chunk number 141: Tutorial_Onemap_reduced_version.Rnw:1328-1329 (eval = FALSE)
###################################################
## (LG3.f2.final<-make.seq(temp.try, 4))


###################################################
### code chunk number 142: Tutorial_Onemap_reduced_version.Rnw:1342-1343 (eval = FALSE)
###################################################
## maps.list<-list(LG1.f2.final, LG2.f2.final, LG3.f2.final)


###################################################
### code chunk number 143: Tutorial_Onemap_reduced_version.Rnw:1348-1349 (eval = FALSE)
###################################################
## draw.map(maps.list, names= TRUE, grid=TRUE, cex.mrk=0.7)


###################################################
### code chunk number 144: Tutorial_Onemap_reduced_version.Rnw:1354-1355 (eval = FALSE)
###################################################
## draw.map(LG1.f2.final, names= TRUE, grid=TRUE, cex.mrk=0.7)


###################################################
### code chunk number 145: Tutorial_Onemap_reduced_version.Rnw:1372-1373 (eval = FALSE)
###################################################
## write.map(maps.list, "fake.f2.onemap.map")


###################################################
### code chunk number 146: Tutorial_Onemap_reduced_version.Rnw:1380-1381 (eval = FALSE)
###################################################
## install.packages("qtl")


###################################################
### code chunk number 147: Tutorial_Onemap_reduced_version.Rnw:1386-1387 (eval = FALSE)
###################################################
## library("qtl")


###################################################
### code chunk number 148: Tutorial_Onemap_reduced_version.Rnw:1392-1394 (eval = FALSE)
###################################################
## raw.file<-paste(system.file("example",package="onemap"),
##                 "fake.f2.onemap.raw", sep="/")


###################################################
### code chunk number 149: Tutorial_Onemap_reduced_version.Rnw:1398-1399 (eval = FALSE)
###################################################
## fake.f2.qtl <- read.cross("mm", file=raw.file, mapfile="fake.f2.onemap.map")


###################################################
### code chunk number 150: Tutorial_Onemap_reduced_version.Rnw:1405-1406 (eval = FALSE)
###################################################
## newmap <- est.map(fake.f2.qtl, tol=1e-6, map.function="kosambi")


###################################################
### code chunk number 151: Tutorial_Onemap_reduced_version.Rnw:1411-1412 (eval = FALSE)
###################################################
## plot.map(fake.f2.qtl, newmap)


###################################################
### code chunk number 152: Tutorial_Onemap_reduced_version.Rnw:1419-1423 (eval = FALSE)
###################################################
## fake.f2.qtl <- calc.genoprob(fake.f2.qtl, step=2)
## out.em <- scanone(fake.f2.qtl, method="em")
## out.hk <- scanone(fake.f2.qtl, method="hk")
## plot(out.em, out.hk, col=c("blue","red"))


###################################################
### code chunk number 153: Tutorial_Onemap_reduced_version.Rnw:1430-1431 (eval = FALSE)
###################################################
## write.cross(fake.f2.qtl, format="qtlcart", filestem="fake.f2.onemap")


###################################################
### code chunk number 154: Tutorial_Onemap_reduced_version.Rnw:1522-1523 (eval = FALSE)
###################################################
## def.rf.3pts(example, "M18", "M8", "M13")


###################################################
### code chunk number 155: Tutorial_Onemap_reduced_version.Rnw:1530-1532 (eval = FALSE)
###################################################
## def.rf.3pts(example, "M18", "M8", "M13", LOD=10, max.rf=0.4)
## def.rf.3pts(example, "M18", "M8", "M13", max.rf=0.4, max.nolink=0.60)


###################################################
### code chunk number 156: Tutorial_Onemap_reduced_version.Rnw:1537-1540 (eval = FALSE)
###################################################
## def.rf.3pts(example, "M18", "M8", "M13")
## def.rf.3pts(example, "M8", "M13", "M7")
## def.rf.3pts(example, "M13", "M7", "M22")


