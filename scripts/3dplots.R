## jobid <- as.numeric(commandArgs(trailingOnly=TRUE))
## if(length(jobid)==0){jobid <- 1}

#pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
#library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
#library('foreach')
#library('LaplacesDemon')
#library('RNetCDF')
#library('Rmpfr')
library('rmarkdown')
#options(bitmapType='cairo')
## mypurpleblue <- '#4477AA'
## myblue <- '#66CCEE'
## mygreen <- '#228833'
## myyellow <- '#CCBB44'
## myred <- '#EE6677'
## myredpurple <- '#AA3377'
## mygrey <- '#BBBBBB'
## mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
## palette(mycolours)
## barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
## barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
## dev.off()

subsamplesizecov <- 2^11
subsamplesizespk <- 2^11
print(paste0('subsample size coverage: ',subsamplesizecov))
print(paste0('subsample size spikes: ',subsamplesizespk))

animaldirs <- list.files(pattern='OT-[0-9]*',include.dirs=TRUE)

for(diritem in animaldirs){
    print(diritem)
    animalid <- unlist(strsplit(diritem,'-'))[2]

    savedir <- paste0('3dplots_ss_',subsamplesizecov,'_',subsamplesizespk,'_',animalid,'/')
    dir.create(savedir)

    sessionfiles <- list.files(path=diritem,pattern='[0-9]*[.]csv')

    for(fileitem in sessionfiles){
        #print(fileitem)
        sessionid <- unlist(strsplit(fileitem,'[.]'))[1]

        neuronfiles <- list.files(path='data_nc/',pattern=paste0('neu_[0-9]*_',animalid,'-',sessionid))
        #print(length(neuronfiles))

        for(neuronitem in neuronfiles){
            ##print(neuronitem)
            fileparts <- unlist(strsplit(neuronitem,'[._]'))
            neuronid <- fileparts[2]
            sessiontype <- fileparts[4]

            outputfile <- paste0(savedir,neuronid,'_',animalid,'_',sessionid,'_',sessiontype,'.htm')

            rmarkdown::render('3dgenerator_XYH.Rmd', output_file=outputfile,
                              params=list(animalid=animalid,
                                          neuronid=neuronid,
                                          sessionid=sessionid,
                                          sessiontype=sessiontype,
                                          subsamplesizecov=subsamplesizecov,
                                          subsamplesizespk=subsamplesizespk
                                          ))
        }
    }
}

print('Done')
