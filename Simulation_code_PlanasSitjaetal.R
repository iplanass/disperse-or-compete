#last modification: 2022/05/11 - Tested on R 3.5 - R 4.1
#This version uses previously generated rasters, instead of generating them itself. 
#Using library(NLMR) could generate a different raster at the beginning of each simulation.
library(ggplot2)
library(stats)
library(sp)
library(RColorBrewer)
require(gridExtra)
library(reshape2)
library(doMC)
library(gridGraphics)
#~ library(NLMR)
library(raster)
library(lattice)
library(landscapetools)
options(error=function()traceback(2))
ptm0<-proc.time()

#******************************#
#------------------------------#
#--- MODIFY PARAMETERS HERE ---#
#------------------------------#
#******************************#

#-----------------#
# NUMBER OF CORES #
#-----------------#

options(cores=8)
registerDoMC()

################################################
############ INITIAL VALUES ####################
################################################

#SIMULATION PARAMETERS
repetitions=8
time=100 # time of simulation
matrix_size=30 # matrix size
t_dyn=10 # Every how many years do we keep values for dynamics
dyn_breaks=10 #breaks for histogram dynamics

#LANDSCAPE PARAMETERS
#~ K=1000 # maximum quality (use for single simulation, instead of K_list below)
SpC=1.2 # Spatial correlation index; if random SpC == 0.001
SEP=0.05 # Stochastic Extinction events Probability
SEPLocal=1 # if 1=local extinction; if 0=continuous disturbance of season
Heterogeneity=1 # 1 = K varies among patches; 0 = Same K for every patch
Binary_land=1
N_patch=2
rasterrandom=system("ls /.../Rasters/Random/*",T) #add directory with Rasters random landscape
rasteragg=system("ls /.../Rasters/Aggregated/*",T) #add directory with Rasters aggregated landscape

#STARTING COLONIES PARAMETERS
#~ SC="LL" #Can choose single SC instead of a SC_list below
R=1.2 # R in logistic equation
Fixed_params=1 # 1= Values of So and E are fixed a priori; 0 = random values
N=matrix_size*matrix_size # initial number of colonies, can have different attributes; define E and So
Qsize=20 #initial size

#-----------------#
#      Evolve     #
#-----------------#

MuE=0.001 # probability of mutation E
er=0.1 # Evolutionary rate --> sigma E
MuSo=0.001 # probability of mutation So
sor=0.1 #evolutionary rate --> sigma So

#---------------------------------------------#
#    Equations : Dispersion and Mortality     #
#---------------------------------------------#

# DISPERSION MODE 01 - STRONG TRADE-OFF
#~ alpha=1 #Minimum distance for equation
#~ delta=20 #Maximum distance for equation
#~ d_th=21 # Threshold dispersion
#~ gamma=50 #steepness

# DISPERSION MODE 02 - WEAK TRADE-OFF
alpha=10 #Minimum distance for equation
delta=20 #Maximum distance for equation
#d_th=21 # Not needed here
gamma=0.002 #steepness

# MORTALITY
aa=0.05 #Minimum mortality for equation
mm=0.95 #Maximum mortality for equation
Te=23 #quorum threshold
beta=20 #steepness

SC_list=c("HH", "HL", "LL", "M")
if (Heterogeneity==1){
	K_list=c(100,500,1000,2000)
}else{
	K_list=c(75, 375, 750, 1500)
}
######################
## DEFINE FUNCTIONS ##
######################

d<-function(nSo){ #DISPERSAL
#~ 			alpha+(delta-alpha)*(nSo^-gamma)/(d_th^-gamma+nSo^-gamma)#DISPERSAL MODE 01
			alpha+(delta-alpha)*exp(-nSo*gamma) #DISPERSAL MODE 02
		}

m<-function(nSo){ #MORTALITY
			aa+(mm-aa)*(nSo^-beta)/(Te^-beta+nSo^-beta)
		}

#********************************************************#
#--------------------------------------------------------#
# 		END OF PARAMETER SECTION		 #
# this code will store results on your working directory #
#--------------------------------------------------------#
#********************************************************#

StartCond<-function(SC){ #DEFINE EI SOI DEPENDING ON INITIAL CONDITION
	if(SC=="HH"){
		EI=c(rep(as.integer(K-5),N)) # Initial energy investment value
		SoI=c(rep(as.integer(K-50),N))# Initial queen propagule size
	}else if(SC=="HL"){
		EI=c(rep(as.integer(K-5),N)) # Initial energy investment value
		SoI=c(rep(5,N))# Initial queen propagule size
	}else if(SC=="LL"){
		EI=c(rep(50,N)) # Initial energy investment value
		SoI=c(rep(5,N))# Initial queen propagule size
	}else if(SC=="M"){
		EI=c(rep(as.integer(K/2),N)) # Initial energy investment value
		SoI=c(rep(as.integer(K/2 - 45),N))# Initial queen propagule size
	}else{
		abort("Starting condition cannot be computed. Please introduce HH, HL, LL or M")
	}
	if(any(SoI<1 | EI<1)==T){
		abort("K value is too low. Please use a higher K.")
	}
	resSC=list(e=EI, s=SoI)
	return(resSC)
}

for (K in K_list){
	for (SC in SC_list){
		if (K==75 & SC=="M"){
			next
		}
		else if (K==100 & SC=="M"){
			next
		}
		
		EI=StartCond(SC)$e #return only e from resSC
		SoI=StartCond(SC)$s #return only s from resSC

		#FOLDERS
		realdate=system("date +%F",T)
		if (Heterogeneity==1){
			if (SpC==0.001){
				newfoldername=paste(realdate,"_R_K",K,"_E",EI[1],"_So",SoI[1],"_rep",repetitions,sep="")
			}else{
				newfoldername=paste(realdate,"_A_K",K,"_E",EI[1],"_So",SoI[1],"_rep",repetitions,sep="")
			}
		}else{
			newfoldername=paste(realdate,"_H_K",K,"_E",EI[1],"_So",SoI[1],"_rep",repetitions,sep="")
		}
		
		system(paste("mkdir ", newfoldername))
		system(paste("mkdir ", newfoldername,"/dynamicsE_model",sep=""))
		system(paste("mkdir ", newfoldername,"/dynamicsSo_model",sep=""))
		system(paste("mkdir ", newfoldername,"/results_model",sep=""))

		###############################################
		#### CREATION OF MATRIXES AND PATCHES #########
		###############################################
		
		mum_type=c(seq(1,N,1))
		columns=c('size','So','E','mum','comp','age','t')
		columns2=c('size','So','E','mum')
		repet=seq(1,repetitions,1)
		ReshE=matrix(0,floor(time/t_dyn)+1,(K/dyn_breaks)+1)
		ReshSo=matrix(0,floor(time/t_dyn)+1,(K/dyn_breaks)+1)
		colsResh=seq(0,(K/dyn_breaks),1)*dyn_breaks
		rowsResh=seq(0,time,t_dyn)
		
		colnames(ReshE)=colsResh
		colnames(ReshSo)=colsResh
		rownames(ReshE)=rowsResh
		rownames(ReshSo)=rowsResh
		
		xy <- expand.grid(seq_len(matrix_size), seq_len(matrix_size))
		
		#####################################################
		###########      SIMULATION      ####################
		#####################################################
		
		RES_TOT<-foreach (rep=1:repetitions, .combine=rbind, .multicombine=TRUE) %dopar%{
			print(rep)
			ptm<-proc.time()
		
			#------------------------------------------------#
			#			 RESOURCES DISTRIBUTION 			 #
			#------------------------------------------------#
			
			if (Heterogeneity==1){
				if (SpC==0.001){
					rland<-read.table(rasterrandom[rep+2],h=T)
				}else{
					rland<-read.table(rasteragg[rep+2],h=T)
				}
				coordinates(rland)<-~ x + y
				gridded(rland)<- TRUE
				mapdist<-raster(rland)
				
				if (Binary_land==1){
					mapdist<-util_classify(x=mapdist,weighting=rep(1/N_patch,N_patch))
				}
				sim1<-as.data.frame(flip(mapdist,2))
				rland=cbind(xy,sim1)
				names(rland) <- c("x","y","K") 
				rland$K=as.integer(rland$K/max(rland$K)*K)
				
		#~ 		#PLOT RASTER GENERATED
		#~ 		coordinates(rland)<-~ x + y
		#~ 		gridded(rland)<- TRUE
		#~ 		rasterDF<-raster(rland)
		#~ 		plot(rasterDF)
		#~ 		#SAVE RASTER
		#~ 		writeRaster(mapdist, filename="multilayer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
		#~ 		mystack = stack("multilayer.tif")
		
			}else{
				
				sim1=rep(K,length(xy[,1]))
				rland=cbind(xy,sim1)
				names(rland) <- c("x","y","K")
			}
				
			#-----------------------#
			#     RESET VALUES      #
			#-----------------------#
			
			t=1 #index for results
			t2=2 #index for size dyn results  
			blobs=matrix(0,length(xy[,1]),length(columns)) # Blobs matrix
			newco=matrix(0,1,length(columns2)) #reproducing blobs
			colnames(blobs)=columns
			colnames(newco)=columns2
			blobs=as.matrix(cbind(rland,blobs))
			newco=as.matrix(cbind(rland[1,],newco))
			
			newco=newco[-1,]
			
			pos=sample.int(length(xy[,1]),N)
			
			E=EI # Initial energy investment value
			So=SoI# Initial queen propagule size
			No=0 #'No' will be E/So
			ReshE[]=0
			ReshSo[]=0
			
			#################
			
			if (Fixed_params==1){
				blobs[pos,'size']=Qsize
				blobs[pos,'E']=E
				blobs[pos,'So']=So
				blobs[pos,'mum']=mum_type
				blobs[pos,c('age','t')]=matrix(c(1,0),N,2,byrow=T)
			}
		
			if (Fixed_params==0){
				for (i in Nvector){
					pos_t=pos[i]
					blobs[pos_t,'size']=Qsize
					blobs[pos_t,'E']=sample(K,1)
					blobs[pos_t,'So']=sample(blobs[pos_t,'E'],1)
					blobs[pos_t,c('age','mum','t')]=c(1,i,0)
				}
				EI=blobs[,"E"]
				SoI=blobs[,"So"]
			}
		
			#SIZE DYN MATRIX TO SAVE DYNAMIC DATA ABOUT SIZE, E & So OF DIFFERENT PATCHES OVER TIME
			if (Heterogeneity==0){
				N_patch=1
				Binary_land=1
			}
			if (Binary_land==0){
				N_patch=1
			}
			
			size_dyn<-matrix(0,(N_patch*7)+3,floor(time/t_dyn)+2)
			colnames(size_dyn)=c("K",rowsResh)
			k_rows=rep(seq(1,N_patch),7)
			k_rows=k_rows[order(k_rows)]
			size_dyn[,"K"]=c(k_rows,0,0,0)
			sd_rows=c(rep(c('E','SDE','So','SDSo','Size','SDSize','Np'),N_patch),'Tsize','SDTsize','mum')
			rownames(size_dyn)=sd_rows
			kp=unique(blobs[,'K']) #generate vector for dyn size results
			kp=kp[order(kp,decreasing=T)]
		
			for (kt in seq(1,N_patch)){ #reset dyn_size results
				size_dyn[which(size_dyn[,"K"]==kt),"K"]=kp[kt]
			}
		#---------------------------------------------------#
		#		SIMULATION	    		    #
		#---------------------------------------------------#
		
			for (i in seq(0,time,1)){
				#-----------------------------------------------------------#
				# 		PHASE 1 - GROWTH, COMPETITION		    #
				#-----------------------------------------------------------#
				
				#Organisms without competition:
				Nvector=which(blobs[,"size"]>0 & blobs[,"comp"]==0)
				if (length(Nvector)>0){
					blobs[Nvector,"size"]=floor(blobs[Nvector,"size"]*exp(R*(1-(blobs[Nvector,"size"]/blobs[Nvector,"K"])))) #logistic equation
				}		
		
				#Competition:
				if ((length(newco[,1])>0)==TRUE){#if any unit is in competition			
					incomp=which(blobs[,"comp"]==1) #resident unit in competition
					newco=rbind(newco,blobs[incomp,c(1:7)])
					newco=newco[order(newco[,"x"],newco[,"y"]),]#order
					
					change=c(1+which(diff(newco[,"x"])!=0 | diff(newco[,"y"])!=0),length(newco[,1]))
					c0=1
					
					for (c in change){
						bch=which(blobs[,"x"]==newco[c0,"x"] & blobs[,"y"]==newco[c0,"y"]) #select unit to replace in blobs matrix
						newco[c0:(c-1),"size"]=floor(newco[c0:(c-1),"size"]*exp(R*(1-(sum(newco[c0:(c-1),"size"])/newco[c0:(c-1),"K"]))))
						alive=which(newco[c0:(c-1),"size"]>1)
						if (length(alive>0)){ #Organisms below size 1 will die and will be eliminated & removed competition
							rprob=newco[c0:(c-1),"size"][alive]/sum(newco[c0:(c-1),"size"][alive])
							survive=which(rmultinom(1,1,rprob)>0)
							blobs[bch,4:7]=newco[(c0+alive[survive]-1),4:7]
							blobs[bch,"age"]=1
							c0=c
						}else{
							c0=c
						}						
					}
					newco=newco[-(1:nrow(newco)),]
				}
				blobs[,"comp"]=0
				blobs[,"t"]=i #register last time step
				
				#-------#
				# SEP	#
				#-------#	
				
				SEPtemp=runif(length(blobs[,1]))-SEP
				blobs[which(SEPtemp<0),"size"]=0
				
				#--------------------------------------#
				# PHASE 2 - Dead organisms are removed #
				#--------------------------------------#
				if (any(blobs[,"size"]<1)==TRUE){ #Organisms below size 1 will die and will be eliminated
					die=which(blobs[,"size"]<1)
					blobs[die,-(1:3)]<-0
				}
						
				Nvector=which(blobs[,"size"]>0)
				Nt=length(Nvector)
				
				if (Nt==0){
					print(paste(rep,"all organisms are dead"))
					hE=rep(0,(K/dyn_breaks)+1)
					hSo=rep(0,(K/dyn_breaks)+1)
					blobsadd<-mapply(rep,c(K,Heterogeneity,SpC,0,0),length(blobs[,1]))
					colnames(blobsadd)=c("Ki","H","SpC","Ei","Soi")
					blobsadd[1:N,c("Ei","Soi")]<-c(EI,SoI)
					result=cbind(blobs,blobsadd)
					write.table(result,file=paste("./",newfoldername,"/results_model/blobs",rep,"results.csv",sep=""),row.names=T,col.names=NA)
					write.table(size_dyn,file=paste("./",newfoldername,"/results_model/size_dyn",rep,".csv",sep=""),row.names=T,col.names=NA)
					write.table(ReshE,file=paste("./",newfoldername,"/dynamicsE_model/dyn",rep,"results.csv",sep=""),row.names=T,col.names=NA)
					write.table(ReshSo,file=paste("./",newfoldername,"/dynamicsSo_model/dyn",rep,"results.csv",sep=""),row.names=T,col.names=NA)
					return(list(hE,hSo,ReshE,ReshSo))
					break
				}
				
				if (any(i==seq(0,time,t_dyn))==T){ #dyamic results
					Telap<-proc.time()[3]-ptm[3]
					Tmiss<-(((Telap*time)/i)-Telap)/3600
					Mins=(Tmiss%%1)*60
					print(paste("[",rep,"]",i,"gen;",round(Tmiss),"h :",sprintf("%.0f", Mins),"min"))
					
					if (Binary_land==0){
					}else{
						for (kt in kp){
							size_dyn[which(rownames(size_dyn)=="E" & size_dyn[,"K"]==kt),t2]=mean(subset(blobs[,"E"],blobs[,"E"]>0 & blobs[,"K"]==kt))
							size_dyn[which(rownames(size_dyn)=="SDE" & size_dyn[,"K"]==kt),t2]=sd(subset(blobs[,"E"],blobs[,"E"]>0 & blobs[,"K"]==kt))
							size_dyn[which(rownames(size_dyn)=="So" & size_dyn[,"K"]==kt),t2]=mean(subset(blobs[,"So"],blobs[,"So"]>0 & blobs[,"K"]==kt))
							size_dyn[which(rownames(size_dyn)=="SDSo" & size_dyn[,"K"]==kt),t2]=sd(subset(blobs[,"So"],blobs[,"So"]>0 & blobs[,"K"]==kt))
							size_dyn[which(rownames(size_dyn)=="Size" & size_dyn[,"K"]==kt),t2]=mean(subset(blobs[,"size"],blobs[,"size"]>0 & blobs[,"K"]==kt))
							size_dyn[which(rownames(size_dyn)=="SDSize" & size_dyn[,"K"]==kt),t2]=sd(subset(blobs[,"size"],blobs[,"size"]>0 & blobs[,"K"]==kt))
							size_dyn[which(rownames(size_dyn)=="Np" & size_dyn[,"K"]==kt),t2]=length(subset(blobs[,"size"],blobs[,"size"]>0 & blobs[,"K"]==kt))
						}
						size_dyn[which(rownames(size_dyn)=="Tsize"),t2]=mean(subset(blobs[,"size"],blobs[,"size"]>0))
						size_dyn[which(rownames(size_dyn)=="SDTsize"),t2]=sd(subset(blobs[,"size"],blobs[,"size"]>0))
						size_dyn[which(rownames(size_dyn)=="mum"),t2]=length(unique(subset(blobs[,"mum"],blobs[,"mum"]>0)))
						t2=t2+1
					}
					hE<-hist(subset(blobs[,"E"],blobs[,"E"]<=K & blobs[,"E"]>0),breaks=seq(0,K+dyn_breaks,dyn_breaks),plot=F)
					hSo<-hist(subset(blobs[,"So"],blobs[,"So"]<=K & blobs[,"So"]>0),breaks=seq(0,K+dyn_breaks,dyn_breaks),plot=F)
					ReshE[t,]= hE$counts #results of histogram to save dynamic of each simulation
					ReshSo[t,]= hSo$counts #results of histogram to save dynamic of each simulation
					t=t+1 #increase index column
					
					blobsadd<-mapply(rep,c(K,Heterogeneity,SpC,0,0),length(blobs[,1]))
					colnames(blobsadd)=c("Ki","H","SpC","Ei","Soi")
					blobsadd[1:N,c("Ei","Soi")]<-c(EI,SoI)
					
					result=cbind(blobs,blobsadd)
					write.table(result,file=paste("./",newfoldername,"/results_model/blobs",rep,"results.csv",sep=""),row.names=T,col.names=NA)
					write.table(size_dyn,file=paste("./",newfoldername,"/results_model/size_dyn",rep,".csv",sep=""),row.names=T,col.names=NA)
					write.table(ReshE,file=paste("./",newfoldername,"/dynamicsE_model/dyn",rep,"results.csv",sep=""),row.names=T,col.names=NA)
					write.table(ReshSo,file=paste("./",newfoldername,"/dynamicsSo_model/dyn",rep,"results.csv",sep=""),row.names=T,col.names=NA)
				}
				
				#---------------------------------------#
				# PHASE 3 - REPRODUCTION AND DISPERSION #
				#---------------------------------------#
				if (i<time){
					for (n in Nvector){
						if (blobs[n,"size"]>=blobs[n,"E"] & blobs[n,"So"]<=blobs[n,"E"]){
							nSo=blobs[n,"So"]
							nE=blobs[n,"E"]
							nNo=floor(nE/nSo)# Number of propagules
							nm=1-m(nSo) #mortality of dispersion per propagule
							blobs[n,"size"]=blobs[n,"size"]-(nSo*nNo) #size of Reproductive unit loses size because of reproduction ####
							PtempN=runif(nNo)
							directionT=sample(seq(0,359,1),nNo,replace=T)
							ndistT=rpois(n=nNo,lambda=d(nSo))
							for (r in seq(1,nNo,1)){
								if (PtempN[r]<nm){  # test survival of new propagule
									ndist=ndistT[r] #Distance to travel each propagule
									direction=directionT[r] #direction to travel
									#--- 0 - 90 degrees ---#
									if (direction<90){
										if(direction==0){
											new_colx=blobs[n,"x"]
											new_coly=blobs[n,"y"]+ndist
										}else{
											newx=floor((sin((direction)*pi/180)*ndist))
											new_colx=blobs[n,"x"]+newx
											new_coly=blobs[n,"y"]+floor(sqrt(ndist**2-newx**2))
										}
									#--- 90 - 180 degrees ---#
									}else if (direction>=90 & direction<180){
										if(direction==90){
											new_colx=blobs[n,"x"]+ndist
											new_coly=blobs[n,"y"]
										}else{
											newx=floor((sin((90-(direction-90))*pi/180)*ndist)) #90-... because opposite angle
											new_colx=blobs[n,"x"]+newx
											new_coly=blobs[n,"y"]-floor(sqrt(ndist**2-newx**2))
										}
									#--- 180 - 270 degrees ---#
									}else if (direction>=180 & direction<270){
										if(direction==180){
											new_colx=blobs[n,"x"]
											new_coly=blobs[n,"y"]-ndist
										}else{
											newx=floor((sin((direction-180)*pi/180)*ndist))
											new_colx=blobs[n,"x"]-newx
											new_coly=blobs[n,"y"]-floor(sqrt(ndist**2-newx**2))
										}
									#--- 270 - 360 degrees ---#
									}else if (direction>=270 & direction<360){
										if(direction==270){
											new_colx=blobs[n,"x"]-ndist
											new_coly=blobs[n,"y"]
										}else{
											newx=floor((sin((90-(direction-270))*pi/180)*ndist)) #90- .. because opposite angle
											new_colx=blobs[n,"x"]-newx
											new_coly=blobs[n,"y"]+floor(sqrt(ndist**2-newx**2))
										}
									}
									
									if (new_colx>matrix_size){#toroid
										new_colx=as.integer(((new_colx/matrix_size)-floor(new_colx/matrix_size))*matrix_size)
									}else if(new_colx<1){
										new_colx=as.integer(((new_colx/matrix_size)-floor(new_colx/matrix_size))*matrix_size)
									}
									if(new_colx==0){
										new_colx=matrix_size
									}
									
									if (new_coly>matrix_size){#toroid
										new_coly=as.integer(((new_coly/matrix_size)-floor(new_coly/matrix_size))*matrix_size)
									}else if(new_coly<1){
										new_coly=as.integer(((new_coly/matrix_size)-floor(new_coly/matrix_size))*matrix_size)
									}
									if(new_coly==0){
										new_coly=matrix_size
									}
											#------------------------#
											#       Evolution	 #
											#------------------------#
									
									newpatch=which(blobs[,"y"]==new_coly & blobs[,"x"]==new_colx) #new patch
									newK=blobs[newpatch,"K"]							
									if (SEPLocal==0){
										if (SEPtemp[newpatch]<0){ #if patch is dead because of SEP
											survive=0
										}else{#patch is not dead zone
											if (blobs[newpatch,"size"]>0){ #if patch is occupied
												survive=1 #goes to newco
											}else{
												survive=2 #establishes in patch
											}
										}
									}else{#organisms can establish in disturbed patches
										if (blobs[newpatch,"size"]>0){ #if patch is occupied
											survive=1 #goes to newco
										}else{
											survive=2 #establishes in patch
										}
										Ptemp=runif(2)
										# Parameter E
										if (Ptemp[1]<MuE){#if probability is lower than Mutation probability, E evolves
											newEtemp=round(sample(seq(nE-er*nE,nE+er*nE,er),1)) #according to sigma --> er
											if (newEtemp>nE){ #if sample value is bigger than nE, we upper round
												newE=ceiling(newEtemp)
											}else{# if sample value is smaller than nE, we lower round, to give a chance of evolving to little propogules, else they will remain of E or So = 1
												newE=floor(newEtemp)
											}
											if (newE<1){
												newE=1
											}
										}else{#parameter do not evolve
											newE=nE
										}
										# Parameter So
										if (Ptemp[2]<MuSo){#if probability is lower than Mutation probability, So evolves
											newSotemp=sample(seq(nSo-sor*nSo,nSo+sor*nSo,sor),1)#according to sigma --> sor
											if (newSotemp>nSo){#if sample value is bigger than nSo, we upper round
												newSo=ceiling(newSotemp)
											}else{# if sample value is smaller than nSo, we lower round, to give a chance of evolving to little propogules, else they will remain of E or So = 1
												newSo=floor(newSotemp)
											}
											if (newSo<1){
												newSo=1
											}
										}else{
											newSo=nSo
										}
										if (survive==2){ #patch was not occupied
											blobs[newpatch,4:7]=c(nSo,newSo,newE,blobs[n,"mum"])
										}else if (survive==1){ #patch occupied; go to newco and resident organisms competition
											newco=rbind(newco,c(new_colx,new_coly,newK,nSo,newSo,newE,blobs[n,"mum"]))
											blobs[newpatch,"comp"]=1
										}
										survive=0
									}
								}else{ #Ptemp>nm ; propagule don't survive
								}
							} #end of reproduction
						}else{#IF PROPAGULE SIZE NOT BIG ENOUGH
						}
					} 
					#compute new lenght of matrix blobs because new organisms were added; age=age+1 for all organisms
					Nvector=which(blobs[,"size"]>0)
					blobs[Nvector,"age"]=blobs[Nvector,"age"]+1
				}
			}
		
			hE<-hist(subset(blobs[,"E"],blobs[,"E"]<=K & blobs[,"E"]>0),breaks=seq(0,K+dyn_breaks,dyn_breaks),plot=F)
			hSo<-hist(subset(blobs[,"So"],blobs[,"So"]<=K & blobs[,"So"]>0),breaks=seq(0,K+dyn_breaks,dyn_breaks),plot=F)
		
			#save files
			blobsadd<-mapply(rep,c(K,Heterogeneity,SpC,0,0),length(blobs[,1]))
			colnames(blobsadd)=c("Ki","H","SpC","Ei","Soi")
			blobsadd[1:N,c("Ei","Soi")]<-c(EI,SoI)
			
			result=cbind(blobs,blobsadd)
			write.table(result,file=paste("./",newfoldername,"/results_model/blobs",rep,"results.csv",sep=""),row.names=T,col.names=NA)
			write.table(size_dyn,file=paste("./",newfoldername,"/results_model/size_dyn",rep,".csv",sep=""),row.names=T,col.names=NA)
			write.table(ReshE,file=paste("./",newfoldername,"/dynamicsE_model/dyn",rep,"results.csv",sep=""),row.names=T,col.names=NA)
			write.table(ReshSo,file=paste("./",newfoldername,"/dynamicsSo_model/dyn",rep,"results.csv",sep=""),row.names=T,col.names=NA)
			print(paste(rep,sprintf("%.2f",proc.time()[3]-ptm[3]))) #print rep number & time with only 2 decimals
			return(list(hE$counts,hSo$counts,ReshE,ReshSo))
		}

		ReshEend=do.call(rbind,RES_TOT[,1])
		ReshSoend=do.call(rbind,RES_TOT[,2])
		write.table(ReshEend,file=paste("./",newfoldername,"/dynamicsE_model/dyn_Hist_results.csv",sep=""),row.names=T,col.names=NA)
		write.table(ReshSoend,file=paste("./",newfoldername,"/dynamicsSo_model/dyn_Hist_results.csv",sep=""),row.names=T,col.names=NA)
		
		#------------------------------------------------#
		#1- Hist of E and So at the end of simulations   #
		#------------------------------------------------#
		
		ResEsum=colSums(ReshEend)
		ResSosum=colSums(ReshSoend)
		
		labels=seq(1,length(ResEsum),1)*dyn_breaks
		fig1=data.frame(gens=labels, valuesE=ResEsum, valuesSo=ResSosum)
		f1<-ggplot(fig1, aes(x=gens, y=valuesE)) + geom_bar(stat = "identity") + labs(x="E", y="Observations")
		f2<-ggplot(fig1, aes(x=gens, y=valuesSo)) + geom_bar(stat = "identity") + labs(x="So", y="")
		
		
		#----------------#
		#2- Dyamics plot #
		#----------------#
		
		ReshEtot=Reduce('+',RES_TOT[,3])/repetitions
		ReshSotot=Reduce('+',RES_TOT[,4])/repetitions
		ReshEtot2=melt(ReshEtot)
		ReshSotot2=melt(ReshSotot)
		colnames(ReshEtot2)=c("generations","E","value")
		colnames(ReshSotot2)=c("generations","So","value")
		f3<-ggplot(ReshEtot2, aes(x = generations, y = E)) + geom_tile(aes(fill = value))+scale_fill_gradient(low="white",high="black")
		f4<-ggplot(ReshSotot2, aes(x = generations, y = So)) + geom_tile(aes(fill = value))+scale_fill_gradient(low="white",high="black")
		
		pdf ("summary.pdf", 16, 8)
		grid.arrange(f1, f2, f3, f4, nrow = 2, top=paste(format(Sys.time(), "%Y %m %d"),", Mu=",MuE,", K=",K,", Ei=",EI[1],", So=",SoI[1],"SpC=",SpC,", gens=",time))
		dev.off()
		
		system(paste("mv summary.pdf",newfoldername))
		print(proc.time()-ptm0)
	}
}
