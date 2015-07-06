# Libraries
library(Rniftilib)
library(plyr)
library(XML)

# Helper Functions
vcat_fun <- function(verbose) {
    if (verbose) {
      fun <- function(msg, ..., newline=TRUE) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
      }
    } else {
      fun <- function(x) invisible(NULL)
    }
    return(fun)
}

read_nifti_as_vector <- function(fname, mask=NULL) {
  nim <- nifti.image.read(fname)
  out <- as.vector(nim[,,])
  if (!is.null(mask)) out <- out[mask]
  out
}

get.coords <- function(xdim, mask=NULL) {
  indices <- lapply(xdim[1:3], seq) # list indices in each dimension
  coords <- expand.grid(indices)  # rows=voxels, cols=x,y,z
  colnames(coords) <- c("i", "j", "k")
  
  if (!is.null(mask)) {
    coords <- coords[mask,]
    rownames(coords) <- 1:nrow(coords)
  }
  
  return(coords)
}

ijk2xyz <- function(nim, crds.ijk) {
  crds.ijk <- cbind(as.matrix(crds.ijk)-1, 1) # -1 converts to 0-based index
  crds.xyz <- tcrossprod(crds.ijk, nim$qto_xyz)
  crds.xyz <- crds.xyz[,-4]
  colnames(crds.xyz) <- c("x", "y", "z")
  as.data.frame(crds.xyz)
}

capitalize <- function(s) {
  paste(toupper(substring(s,1,1)), substring(s,2),
        sep="")
}

decapitalize <- function(s) {
  paste(tolower(substring(s,1,1)), substring(s,2),
        sep="")
}

run <- function(raw.cmd, ...) {
  cmd <- sprintf(raw.cmd, ...)
  system(cmd)
}


# Peak Detection
## note: separation distance is in voxels
detect.peaks <- function(statfile, maskfile, peakfile, 
                         smooth.signif=T, sep.dist.voxs=6) 
{
  nim      <- nifti.image.read(statfile)  
  fwhm     <- mean(nim$pixdim) # will smooth by one voxel
  sep.dist <- sep.dist.voxs * mean(nim$pixdim)
  
  # Temporary files
  clustfile  <- tempfile(pattern="clust", fileext=".nii.gz")
  smoothfile <- tempfile(pattern="smooth", fileext=".nii.gz")
  tablefile  <- tempfile(pattern="table", fileext=".txt")
  
  # Smooth stats
  if (smooth.signif) {    
    # Get clusters
    run("3dclust -savemask %s -dxyz=1 1.75 0 %s", clustfile, statfile)
    # Smoothin within significant clusters
    run("3dBlurInMask -input %s -FWHM %.4f -Mmask %s -prefix %s", statfile, fwhm, clustfile, smoothfile)
  } else {
    # Smooth within mask
    run("3dBlurInMask -input %s -FWHM %.4f -mask %s -prefix %s", statfile, fwhm, maskfile, smoothfile)    
  }
  
  # Get the peaks
  run("3dExtrema -prefix %s -maxima -volume -closure -sep_dist %.4f -mask_file %s %s > %s", peakfile, sep.dist, maskfile, smoothfile, tablefile)
  
  # Temporary files
  file.remove(clustfile, smoothfile, tablefile)
  
  return(TRUE)
}


# Read in ROIs
read.curv <- function(mask=NULL) {
  # Read in data
  nim  <- nifti.image.read("CurvVol_2mm.nii.gz")
  img  <- as.vector(nim[,,])
  
  # Make data binary
  img  <- (img>0)*2 + (img<0)*1 # 1 = gyrus, 2 = sulcus
  
  # Create labels
  df   <- data.frame(index=1:2, name=c("gyrus", "sulcus"))
  
  # Mask
  if (!is.null(mask)) img <- img[mask]
  
  # Return
  list(labels=df, data=img)
}
read.brodmann <- function(mask=NULL) {
  # Read in data
  nim  <- nifti.image.read("brodmann_2mm.nii.gz")
  img  <- as.vector(nim[,,])
  
  # Create labels
  labs <- sort(unique(img[img!=0]))
  df   <- data.frame(index=labs, name=as.character(labs))
  
  # Mask
  if (!is.null(mask)) img <- img[mask]
  
  # Return
  list(labels=df, data=img)
}
read.freesurfer <- function(mask=NULL) {
  # Read Data
  nim <- nifti.image.read("aparc+aseg_2mm.nii.gz")
  img <- as.vector(nim[,,]) # convert data to vector
  
  # Read in the labels
  labs <- read.table("aparc+aseg_final.labels", header=T)
  colnames(labs) <- c("index", "hemi", "name")
  labs$name <- gsub("-", " ", labs$name)
  
  # Mask
  if (!is.null(mask)) img <- img[mask]

  list(labels=labs, data=img)
}
read.cerebellum <- function(mask=NULL) {
  # Read in labels
  data     <- xmlParse("Cerebellum_MNIfnirt.xml")
  xml_data <- xmlToList(data)
  df       <- ldply(xml_data$data, function(x) {
    c(index=as.integer(x$.attrs[1])+1, 
      name=sub("Left |Right ", "", as.character(x$text)), 
      full.name=as.character(x$text))
  })[,-1]
  
  # Read in the data
  img <- as.vector(nifti.image.read("Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz")[,,])
  
  # Mask
  if (!is.null(mask)) img <- img[mask]
  
  # Return
  list(labels=df, data=img)
}
read.yeo <- function(mask=NULL) {
  labels <- c("Visual", "Somatomotor", "Dorsal Attention", "Ventral Attention", "Limbic", "Frontoparietal", "Default")
  df  <- data.frame(
    index=1:length(labels), 
    name=labels
  )
  img <- as.vector(nifti.image.read("7networks_mni152_2mm.nii.gz")[,,])
  
  # Mask
  if (!is.null(mask)) img <- img[mask]

  list(labels=df, data=img)
}

# Parse the peaks to region labels
parse.closest.rois <- function(rois, nim, mask, search.range=5, progress="text", parallel=F) {
    # TODO: check that the dims are the same in the data!
    
    # Get peaks
    peaks <- as.vector(nim[,,])[mask]
    
    # Get coordinates to search
    crds <- get.coords(dim(nim), mask)
    crds.vals <- crds[rois$data!=0,] # can restrict search to areas that exist
    
    # Loop through peaks to get labels
    peak.inds <- which(peaks==1)
    rr <- search.range # voxels
    peak.labs <- ldply(peak.inds, function(peak.ind) {
      if (rois$data[peak.ind]==0) {
        loc <- crds[peak.ind,] # current location
        # Restrict search to a box within search range or rr
        ilim <- crds.vals$i < (rr + loc$i) & crds.vals$i > (rr - loc$i)
        jlim <- crds.vals$j < (rr + loc$j) & crds.vals$j > (rr - loc$j)
        zlim <- crds.vals$k < (rr + loc$k) & crds.vals$k > (rr - loc$k)
        lim  <- ilim & jlim & zlim
        crds.vals.lim <- crds.vals[lim,]
        # Check if anything is within search range
        if (sum(lim)==0) {
          ind <- peak.ind
        } else {
          # Get closest value
          tmp <- sweep(crds.vals.lim, 2, as.numeric(loc))
          d2 <- rowSums(tmp^2) # euclidean squared distance
          closest.ind <- as.numeric(names(which.min(d2))) # the names have the original index
          closest.d <- sqrt(min(d2))
          # Return closest index if within range
          if (closest.d <= rr) ind <- closest.ind
          else ind <- peak.ind
        }
      } else {
        ind <- peak.ind
        closest.d <- 0
      }
      
      val <- rois$data[ind]
      if (val == 0) {
        ret <- c(name="", dist=999)
      } else {
        ret <- c(name=as.character(rois$labels$name[rois$labels$index==val]), dist=closest.d)
      }
      
      ret
    }, .progress=progress, .parallel=parallel)
  
    # Return
    peak.labs
}

read_nifti_as_clust <- function(statfile, mask) {
  clustfile <- tempfile(pattern="clust", fileext=".nii.gz")
  # note 1.75 radius is basically a box around the voxel
  run("3dclust -savemask %s -dxyz=1 1.75 0 %s", clustfile, statfile)
  clust <- read_nifti_as_vector(clustfile, mask)
  file.remove(clustfile)
  return(clust)
}

# Wrapper script
detect.peaks.and.whereami <- function(statfile, maskfile, sep.dist.voxs=6, ...) 
{
  peakfile  <- tempfile(pattern="peaks", fileext=".nii.gz")
  detect.peaks(statfile, maskfile, peakfile, 
               smooth.signif=T, sep.dist.voxs=sep.dist.voxs)
  df.tab <- whereami(statfile, maskfile, peakfile, ...)
  file.remove(peakfile)
  df.tab
}

whereami <- function(statfile, maskfile, peakfile, search.range=5, 
                     progress="text", parallel=F, verbose=T, whereami.dir=NULL) 
{
  vcat <- vcat_fun(verbose)
  
  curdir <- getwd()
  if (is.null(whereami.dir)) {
    whereami.dir <- Sys.getenv("WHEREAMI_DIR")
  }
  if (!file.exists(whereami.dir)) {
    stop("Please set your WHEREAMI_DIR environmental variable or specify as an argument") 
  }
  
  # Read in data
  vcat("reading in data")
  nim.peak <- nifti.image.read(peakfile)
  mask     <- as.logical(read_nifti_as_vector(maskfile))
  peaks    <- read_nifti_as_vector(peakfile, mask)
  stats    <- read_nifti_as_vector(statfile, mask)
  clust    <- read_nifti_as_clust(statfile, mask)
  
  # Get coordinates in xyz (mm space)
  vcat("coordinates in mm")
  crds.ijk <- get.coords(dim(nim.peak), mask)
  crds     <- ijk2xyz(nim.peak, crds.ijk)
  crds.peaks <- crds[peaks==1,]
  
  # Get if left or right. Negative x means left
  vcat("hemisphere")
  hemis    <- factor(sign(crds$x[peaks==1]), levels=c(-1,0,1), labels=c("L", "B", "R"))
  
  # Read in ROIs and parse region associated with each peak
  roi.sets <- c("curv", "freesurfer", "brodmann", "yeo", "cerebellum")
  setwd(whereami.dir)
  df.tab2  <- laply(roi.sets, function(roi.set) {
    vcat("determining %s regions", roi.set)
    lst.rois <- do.call(sprintf("read.%s", roi.set), list(mask))
    df.labs  <- parse.closest.rois(lst.rois, nim.peak, mask, search.range, 
                                   progress=progress, parallel=parallel)
    ret      <- as.character(df.labs[,1])
    ret
  })
  setwd(curdir)
  df.tab2  <- as.data.frame(t(df.tab2))
  colnames(df.tab2) <- roi.sets
  ## only keep cerebellum if freesurfer says cerebellum
  df.tab2$cerebellum[df.tab2$freesurfer!="Cerebellum"] <- ""
  ## lowercase crus and vermis
  df.tab2$cerebellum <- sub("Crus", "crus", df.tab2$cerebellum)
  df.tab2$cerebellum <- sub("Vermis", "vermis", df.tab2$cerebellum)
  ## only keep brodmann area if freesurfer not cerebellum, brain-stem, or subcortical
  ## only keep gyrus/sulcus if freesurfer not cerebellum, brain-stem, or subcortical
  inds <- df.tab2$freesurfer %in% c("Brain Stem", "Cerebellum", "Thalamus", "Caudate", "Putamen", "Pallidum", "Hippocampus", "Amygdala", "Accumbens", "VentralDC")
  df.tab2$brodmann[inds] <- ""
  df.tab2$curv[inds] <- ""
  ## remove sulcus/gyrus from insula and precuneus...TODO
  
  ## remove yeo network in certain subcortical regions
  inds <- df.tab2$freesurfer %in% c("Brain Stem", "Thalamus", "Hippocampus", "Amygdala", "VentralDC")
  df.tab2$yeo[inds] <- ""
  
  # Combine it all together
  vcat("combining it all together")
  region <- gsub(" $", "", paste(capitalize(tolower(df.tab2$freesurfer)), df.tab2$curv, sep=" "))
  region <- gsub("Ventraldc", "VentralDC", region)
  region <- gsub(", $", "", paste(region, df.tab2$cerebellum, sep=", "))
  df.tab <- data.frame(
    Cluster = clust[peaks==1], 
    Region = region,  
    H = hemis, 
    BA = df.tab2$brodmann, 
    x = crds.peaks$x, 
    y = crds.peaks$y, 
    z = crds.peaks$z, 
    Z.score = stats[peaks==1], 
    Network = capitalize(tolower(as.character(df.tab2$yeo)))
  )
  
  # Sort by cluster, y-axis, x-axis, z-axis
  vcat("sort")
  inds <- with(df.tab, order(Cluster, Region, -1*y, x, z))
  df.tab <- df.tab[inds,]
  
  df.tab
}

test.whereami <- function() {
  source("whereami.R")
  statfile <- "test_zstat.nii.gz"
  maskfile <- "MNI152_T1_2mm_brain_mask_dil.nii.gz"
  peakfile <- "test_peaks2.nii.gz"
  file.remove(peakfile)
  detect.peaks(statfile, maskfile, peakfile)
  df.tab <- whereami(statfile, maskfile, peakfile, search.range=5, whereami.dir=getwd())
  print(df.tab)
  file.remove(peakfile)
}
