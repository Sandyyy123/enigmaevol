options(stringsAsFactors=FALSE);
library(plotly);
library(gplots);
#setwd('~/Documents/papers/ENIGMA-Evolution/RegionalPlottingFreesurfer');

##Load in all the Freesurfer surfaces
load('FreesurferRegionalObjs.Rdata')

##Correlations with SDS scores at different thresholds
#Set colorbar values
minZ = -0.0165;
maxZ = 0.0165;
    
fcorvals = (paste0("../SDS_ancreg1KGPphase3_noGC_MA6/SDS_bjk_ancreg_1kblocks.csv"));
corvals = read.csv(fcorvals);
thisSA = corvals[c(grep("surfavg",corvals$X),grep("Mean_Full_SurfArea",corvals$X)),];
thisSA$region = sapply(thisSA$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]});
thisTH = corvals[c(grep("thickavg",corvals$X),grep("Mean_Full_Thickness",corvals$X)),];
thisTH$region = sapply(thisTH$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]});

##Take only the left hemisphere obj files (everything in ENIGMA is mean bilateral)
##So no need for both hemispheres
objs = objs[1:35];
names(objs) = sapply(names(objs),function (x) {unlist(strsplit(x,".",fixed=TRUE))[length(unlist(strsplit(x,".",fixed=TRUE)))]});

colorshere = colorRampPalette(c("midnightblue","white","darkred"))(50);
colorshere[25] = "#BABABC";
##Remove axes
ax = list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)

    thisSA$padj = p.adjust(thisSA$BJK_ESTIM_PVAL,"BH");
    thisSA$BJK_ESTIM_AVE[which(thisSA$padj > 0.05)] = 0;
    ##thisSA$clumpedcorestimate[which(thisSA$padj > 0.05)] = 0;
    ##Make sure the order of each corresponds
    inds = match(names(objs),thisSA$region);
    thisSA = thisSA[inds,];
    #maxZ = max(thisSA$clumpedcorestimate,na.rm=TRUE);
    #minZ = min(thisSA$clumpedcorestimate,na.rm=TRUE);
    ##Ensure the color scale is symmetrical from positive to negative
    #indmax = which.max(c(abs(maxZ),abs(minZ)));
    #if (indmax==1) {
    #    minZ = -abs(maxZ);
    #} else {
    #    maxZ = abs(minZ)
    #}
    ##Loop through each region for plotting
    for (i in 1:length(objs)) {
        ftitle = paste0("SA BJK");
        obj = objs[[i]];
        ##If it is the first plot, set up the mesh
        if (i==1) {
            p = plot_ly(type = 'mesh3d',
                        x = obj$shapes[[1]]$positions[1,],
                        y = obj$shapes[[1]]$positions[2,],
                        z = obj$shapes[[1]]$positions[3,],
                        i = obj$shapes[[1]]$indices[1,],
                        j = obj$shapes[[1]]$indices[2,],
                        k = obj$shapes[[1]]$indices[3,],
                        ##for two of the vertices, do the whole range of colors (you can't see the in the plots, but without it I couldn't get the color ranges to work)
                        intensity = rep(thisSA$BJK_ESTIM_AVE[i],length(obj$shapes[[1]]$positions[1,])),
                        opacity=1,
                        cauto=FALSE,
                        cmax=maxZ,
                        cmin=minZ,
                        colors = colorshere,
                        showscale=TRUE)
        } else {
            p = add_mesh(p = p,  
                         x = obj$shapes[[1]]$positions[1,],
                         y = obj$shapes[[1]]$positions[2,],
                         z = obj$shapes[[1]]$positions[3,],
                         i = obj$shapes[[1]]$indices[1,],
                         j = obj$shapes[[1]]$indices[2,],
                         k = obj$shapes[[1]]$indices[3,],
                         intensity = rep(thisSA$BJK_ESTIM_AVE[i],length(obj$shapes[[1]]$positions[1,])),
                         opacity=1,
                         cauto=FALSE,
                         cmax=maxZ,
                         cmin=minZ,
                         showscale=FALSE)
        }
    }
    p = layout(p=p,scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 2.3, y = 0, z = 0))),title=ftitle);
    htmlwidgets::saveWidget(p, file = paste0("BJK_medial_SA_clump1.html"));
    p = layout(p=p,scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = -2, y = 0, z = 0))),title=ftitle);
    htmlwidgets::saveWidget(p, file = paste0("BJK_lateral_SA_clump1.html"));
    
##Make a brain all of one color for the global SA
fcorvals = (paste0("../SDS_ancreg1KGPphase3_noGC_MA6/SDS_bjk_ancreg_1kblocks.csv"));
corvals = read.csv(fcorvals);
thisSA = corvals[c(grep("surfavg",corvals$X),grep("Mean_Full_SurfArea",corvals$X)),];
thisSA$region = sapply(thisSA$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]});
##Take only the left hemisphere obj files (everything in ENIGMA is mean bilateral)
##So no need for both hemispheres
objs = objs[1:35];
names(objs) = sapply(names(objs),function (x) {unlist(strsplit(x,".",fixed=TRUE))[length(unlist(strsplit(x,".",fixed=TRUE)))]});
    
colorshere = colorRampPalette(c("midnightblue","white","darkred"))(50);
colorshere[25] = "#BABABC";
##Remove axes
ax = list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
)
    
thisSA$padj = p.adjust(thisSA$BJK_ESTIM_PVAL,"BH");
thisSA$BJK_ESTIM_AVE[which(thisSA$padj > 0.05)] = 0;

valueforwholebrain = thisSA[which(thisSA$region=="Full"),]$BJK_ESTIM_AVE;

    ##Loop through each region for plotting
    for (i in 1:length(objs)) {
        ftitle = paste0("Global SA");
        obj = objs[[i]];
        ##If it is the first plot, set up the mesh
        if (i==1) {
            p = plot_ly(type = 'mesh3d',
                        x = obj$shapes[[1]]$positions[1,],
                        y = obj$shapes[[1]]$positions[2,],
                        z = obj$shapes[[1]]$positions[3,],
                        i = obj$shapes[[1]]$indices[1,],
                        j = obj$shapes[[1]]$indices[2,],
                        k = obj$shapes[[1]]$indices[3,],
                        ##for two of the vertices, do the whole range of colors (you can't see the in the plots, but without it I couldn't get the color ranges to work)
                        intensity = rep(valueforwholebrain,length(obj$shapes[[1]]$positions[1,])),
                        opacity=1,
                        cauto=FALSE,
                        cmax=maxZ,
                        cmin=minZ,
                        colors = colorshere,
                        showscale=TRUE)
        } else {
            p = add_mesh(p = p,  
                         x = obj$shapes[[1]]$positions[1,],
                         y = obj$shapes[[1]]$positions[2,],
                         z = obj$shapes[[1]]$positions[3,],
                         i = obj$shapes[[1]]$indices[1,],
                         j = obj$shapes[[1]]$indices[2,],
                         k = obj$shapes[[1]]$indices[3,],
                         intensity = rep(valueforwholebrain,length(obj$shapes[[1]]$positions[1,])),
                         opacity=1,
                         cauto=FALSE,
                         cmax=maxZ,
                         cmin=minZ,
                         showscale=FALSE)
        }
    }
p = layout(p=p,scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = 2.3, y = 0, z = 0))),title=ftitle)
htmlwidgets::saveWidget(p, file = paste0("BJK_medial_GlobalSA_clump1.html"))
p = layout(p=p,scene = list(xaxis=ax,yaxis=ax,zaxis=ax,camera = list(eye = list(x = -2, y = 0, z = 0))),title=ftitle)
htmlwidgets::saveWidget(p, file = paste0("BJK_lateral_GlobalSA_clump1.html"))