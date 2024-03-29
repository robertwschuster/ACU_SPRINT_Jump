#-----------------------------------------------------------------------------------------
#
# Load and analyze squat jump data (*.csv files)
# Robert Schuster (ACU SPRINT)
# August 2022
#
#-----------------------------------------------------------------------------------------


## clear environment
rm(list = ls())


# Functions ------------------------------------------------------------------------------
# https://github.com/stas-g/findPeaks
# a 'peak' is defined as a local maxima with m points either side of it being smaller than it. 
# hence, the bigger the parameter m, the more stringent is the peak finding procedure
find_peaks <- function(x, m = 3) {
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

# moving standard deviation
movstd <- function(vec, width)      
  return(c(rep(NA,width-1), sapply(seq_along(vec[width:length(vec)]),      
                                   function(i) sd(vec[i:(i+width-1)]))))

# find flight as the longest period where Fz is smaller than a threshold
find_flight <- function(fz, th) {
  fz[fz <= th] <- 1
  fz[fz > th] <- 0
  
  z <- rle(fz == 1)
  zmi <- which(z$lengths == max(z$lengths[z$values]))
  
  out <- list(start = sum(z$lengths[seq_len(zmi-1)]) + 1,
              end = sum(z$lengths[seq_len(zmi)]))
  return(out)
}


# Load data ------------------------------------------------------------------------------
# select files to analyze
flist <- choose.files(caption = "Choose squat jump files")

# load data and header separately
skip <- list()
data <- list()
cols <- c("Time","Z.Left","Z.Right","Left","Right")
header <- list()
for (f in flist) {
  fn <- basename(f)
  # determine length of header
  skip[[fn]] <- read.csv(f, nrows = 20, header = F, blank.lines.skip = F, 
                         col.names = paste0("V",seq_len(max(count.fields(f, sep = ',')))), fill = T)
  if (any(skip[[fn]] == "Time", na.rm = T)) {
    skip[[fn]] <- which(skip[[fn]] == 'Time')-1
    
    # load data (without header) and keep only Time, Left and Right columns
    data[[fn]] <- read.csv(f, skip = skip[[fn]])
    if (any(names(data[[fn]]) %in% cols, na.rm = T)) {
      data[[fn]] <- data[[fn]][,which(names(data[[fn]]) %in% cols)] # keep only left, right and time columns
      colnames(data[[fn]]) <- gsub('Z.','',colnames(data[[fn]])) # remove 'Z." from column names
      
      data[[fn]]$Total <- data[[fn]]$Left + data[[fn]]$Right
      
      # load header
      header[[fn]] <- read.csv(f, nrows = skip[[fn]]-1, header = F)
      
    } else {
      # print warning but continue execution with other files
      warning(paste("The file ", fn, " does not contain the required data and was not processed"))
      flist <- flist[-which(flist == f)] # remove the file from file list
    }
  } else {
    # print warning but continue execution with other files
    warning(paste("The file", fn, "does not contain the required data and was not processed"))
    flist <- flist[-which(flist == f)] # remove the file from file list
  }
}
rm(skip,cols,f,fn)

fpath <- unique(dirname(flist))
fnames <- basename(flist)

bodymass <- list()
freq <- list()
for (f in fnames) {
  # find weight for force normalisation
  if (any(header[[f]] == "Weight", na.rm = T)) {
    bodymass[[f]] <- as.numeric(na.omit(header[[f]][header[[f]] == 'Weight',2])) * 9.81
  } else {
    bodymass[[f]] <- median(data[[f]]$Total[data[[f]]$Total >= 300])
  }
  # find frequency
  freq[[f]] <- as.numeric(na.omit(header[[f]][header[[f]] == 'Frequency',2]))
  if (freq[[f]] < 1000) {
    message(paste(f,": Sampling frequency =",freq[[f]],"\nFrequencies > 1000 Hz are recommended"))
  }
}
rm(f,header)


# Determine number of repetition per file ------------------------------------------------
fmaxi <- list()
for (f in fnames) {
  # find peaks (at least 2s apart from each other)
  d <- freq[[f]]*2 # frequency * 2s
  pks <- find_peaks(data[[f]]$Total, m = d)
  plot(x = data[[f]]$Time, y = data[[f]]$Total, type = "l", lwd = 2, 
       xlab = "Time [s]", ylab = "Force [N]", main = f)
  points(x = data[[f]]$Time[c(1:nrow(data[[f]]))[pks]], y = data[[f]]$Total[pks], col = "red")
  # determine a cutoff under which all other peaks cannot be considered IMTP trials
  cutoff <- (max(data[[f]]$Total[pks])-bodymass[[f]])/2.75 + bodymass[[f]] # 36% of BM normalised max force
  # cutoff <- max(data[[fn]]$Total[pks])-250 # within 250 N of max value
  abline(h = cutoff, col = "green", lty = "dashed")
  abline(h = bodymass[[f]], col = "blue", lty = "dashed")
  abline(h = 0, col = "red", lty = "dashed")
  # only keep peaks above cutoff and more than 5s apart
  pks <- pks[which(data[[f]]$Total[pks] > cutoff)]
  d <- freq[[f]]*2.5 # frequency * 2.5s
  if (any(diff(pks) < d)) {
    pks <- pks[-(which(diff(pks) < d)+1)]
  }
  points(x = data[[f]]$Time[c(1:nrow(data[[f]]))[pks]], y = data[[f]]$Total[pks], col = "red", pch = 16)
  
  # determine start and end of rep
  for (p in 1:length(pks)) {
    # define arbitrary period to look for start of rep
    if (p > 1) { # if not first rep
      ds <- min(freq[[f]]*3, (pks[p]-pks[p-1])/2) # start of period = halfway between two adjacent peaks or 3S before peak
    } else {
      ds <- min(freq[[f]]*3, pks[p]-1) # else start of period = either start of trial or 3s before peak
    }
    
    # define end of rep
    if (p < length(pks)) { # if not last rep
      de <- min(freq[[f]]*2, (pks[p+1]-pks[p])/2) # end of rep = halfway between two adjacent peaks or 2s after peak
    } else {
      de <- min(freq[[f]]*2, length(data[[f]]$Total)-pks[p]) # else end of rep = either end of trial or 2s after peak
    }
    w = (pks[p]-ds):(pks[p]+de)
    
    if (length(pks) > 1) {
      n <- paste(f,p,sep = '_')
      data[[n]] <- data[[f]][w,]
      
      bodymass[[n]] <- bodymass[[f]]
      freq[[n]] <- freq[[f]]
      
      fmaxi[[n]] <- (pks[p] - w[1]) + 1
      if (p == length(pks)) {
        data[[f]] <- NULL
        fmaxi[[f]] <- NULL
        freq[[f]] <- NULL
        bodymass[[f]] <- NULL
        
      }
    } else {
      data[[f]] <- data[[f]][w,]
      fmaxi[[f]] <- (pks - w[1]) + 1
    }
  }
}
fmaxi <- fmaxi[match(names(data), names(fmaxi))]
rm(f,d,pks,cutoff,p,ds,de,w,n)


# Determine start and end of flight ------------------------------------------------------
th <- numeric(length(data))
sj <- numeric(length(data))
flight <- matrix(0,length(data),2)
colnames(flight) <- c('start','end')
for (f in 1:length(data)) {
  # threshold = 5 SD of 1s weighing period before flight or smallest SD or 1s rolling window
  th[f] <- min(movstd(data[[f]]$Total, (freq[[f]]*1)), na.rm = T)*5
  
  # flight thresholds = 5N start, 20N end
  sft <- 5
  eft <- 20
  
  # start and end of flight = first and last instance when force < flight thresholds
  zero <- find_flight(data[[f]]$Total[1:fmaxi[[f]]], sft)
  flight[f,1] <- min(which(data[[f]]$Total[1:zero$end] < sft))
  flight[f,2] <- min(which(data[[f]]$Total[zero$start:length(data[[f]]$Total)] > eft)) + zero$start
  
  # start of jump movement
  ep <- which(data[[f]]$Total[1:flight[f,1]] == max(data[[f]]$Total[1:flight[f,1]]))
  sj[f] <- max(which(data[[f]]$Total[1:ep] < (bodymass[[f]] + th[f]))) # - (freq[[f]]*0.03)
}
rm(f,sft,eft,ep)


# Check for indicators of poor trial -----------------------------------------------------
ptq <- matrix('_',length(data),2)
colnames(ptq) <- c('Unstable weighing','Countermovement')
rownames(ptq) <- names(data)
for (f in 1:length(data)) {
  sf <- flight[f,1]
  ef <- flight[f,2]
  
  plot(x = data[[f]]$Time, y = data[[f]]$Total, type = "l", lwd = 2, 
       xlab = "Time [s]", ylab = "Force [N]", main = names(data)[f])
  abline(h = bodymass[[f]], col = "green", lwd = 2, lty = "dashed")
  abline(v = data[[f]]$Time[sf], col = "red", lwd = 2, lty = "dotted")
  abline(v = data[[f]]$Time[ef], col = "red", lwd = 2, lty = "dotted")
  abline(v = data[[f]]$Time[sj[f]], col = "blue", lwd = 2, lty = "dotdash")
  
  # unstable weighing period before jump (change in force > 50 N)
  if (any(diff(data[[f]]$Total[1:sf]) > 50)) {
    ptq[f,1] <- 'X'
    # message(paste("Warning:", names(data)[f], "does not have a stable weighing period"))
    # # ask whether to continue processing or not
    # resp <- readline("Would you like to continue processing anyway? (Y/N) ")
    # if (grepl(resp, 'n', ignore.case = TRUE)) {
    #   stop("Processing stopped")
    # }
  }
  
  # countermovement before jump
  if (any(data[[f]]$Total[(sj[f] - freq[[f]]*0.25):sj[f]] < (bodymass[[f]] - th[f]))) {
    ptq[f,2] <- 'X'
    # message(paste("Warning:", names(data)[f], "contains a countermovement prior to the squat jump"))
    # # ask whether to continue processing or not
    # resp <- readline("Would you like to continue processing anyway? (Y/N) ")
    # if (grepl(resp, 'n', ignore.case = TRUE)) {
    #   stop("Processing stopped")
    # }
  }
}
rm(f,sf,ef)


# Extract performance metrics ------------------------------------------------------------
# jump height
FT <-  numeric(length(data))
JH <-  matrix(0,length(data),3)
colnames(JH) <- c('JH_ft','JH_J','JH_di')
for (f in 1:length(data)) {
  # flight time (FT)
  FT[f] <- data[[f]]$Time[flight[f,2]] - data[[f]]$Time[flight[f,1]]
  u <- 9.81 * (FT[f]/2) # velocity at take-off
  JH[f,1] <- (0^2 - u^2) / (2*-9.81)
  # impulse (J)
  Fr <- data[[f]]$Total[sj[f]:flight[f,1]] - bodymass[[f]] # resultant force between start of movement and take-off
  J <- sum(diff(data[[f]]$Time[sj[f]:flight[f,1]]) * (head(Fr,-1) + tail(Fr,-1)))/2 # impulse (trapz integration of resultant force)
  JH[f,2] <- ((J/(bodymass[[f]]/9.81))^2)/(2*9.81)
  # double integration (di)
  a <- (data[[f]]$Total[sj[f]:flight[f,2]] / (bodymass[[f]] / 9.81)) - 9.81 # acceleration between start of movement and take-off
  v <- cumsum(c(0,(a[1:(length(a)-1)] + a[2:length(a)])/2 * diff(data[[f]]$Time[sj[f]:flight[f,2]]))) # velocity (trapz integration of acceleration)
  d <- cumsum(c(0,(v[1:(length(v)-1)] + v[2:length(v)])/2 * diff(data[[f]]$Time[sj[f]:flight[f,2]]))) # displacement (trapz integration of velocity)
  JH[f,3] <- max(d) - d[length(d)]
}
rm(f,u,Fr,J,a,v,d)

# Performance metrics
PM <-  matrix(0,length(data),23)
colnames(PM) <- c('v t-o',
                  'peak force','rel peak force','peak velocity','peak power','rel peak power',
                  'mean force','rel mean force','mean velocity','mean power','rel mean power',
                  'RFD-50','RFD-100','RFD-150','RFD-200','RFD-MM',
                  'J-50','J-100','J-150','J-200','J-T','FT:CT','Contraction T')
for (f in 1:length(data)) {
  # impulse
  fz <- data[[f]]$Total[sj[f]:flight[f,1]] # Fz = force between start of movement and take-off
  j <- cumsum(c(0,(fz[1:(length(fz)-1)] + fz[2:length(fz)])/2 * diff(data[[f]]$Time[sj[f]:flight[f,1]]))) # impulse (trapz integration of force)
  # power
  a <- (data[[f]]$Total[sj[f]:flight[f,2]] / (bodymass[[f]] / 9.81)) - 9.81 # acceleration between start of movement and landing
  v <- cumsum(c(0,(a[1:(length(a)-1)] + a[2:length(a)])/2 * diff(data[[f]]$Time[sj[f]:flight[f,2]]))) # velocity (trapz integration of acceleration)
  p <- data[[f]]$Total[sj[f]:flight[f,2]] * v # power = Fz * velocity
  # peak power index
  ppi <- which(p == max(p))
  # peak & min force indices
  pfi <- which(data[[f]]$Total[sj[f]:flight[f,2]] == max(data[[f]]$Total[sj[f]:flight[f,2]])) + sj[f]
  mfi <- which(data[[f]]$Total[sj[f]:pfi] == min(data[[f]]$Total[sj[f]:pfi])) + sj[f]
  
  PM[f,1] <- v[flight[f,1]-sj[f]] # take-off velocity
  
  PM[f,2] <- max(data[[f]]$Total[sj[f]:flight[f,2]]) # peak force
  PM[f,3] <- max(data[[f]]$Total[sj[f]:flight[f,2]]) / bodymass[[f]] # relative peak force
  PM[f,4] <- max(v) # peak velocity
  PM[f,5] <- max(p) # peak power
  PM[f,6] <- max(p) / bodymass[[f]] # relative peak power
  
  PM[f,7] <- mean(data[[f]]$Total[sj[f]:flight[f,1]]) # mean force
  PM[f,8] <- mean(data[[f]]$Total[sj[f]:flight[f,1]]) / bodymass[[f]] # relative mean force
  PM[f,9] <- mean(v[1:which(v == max(v))]) # mean velocity
  PM[f,10] <- mean(p[1:ppi]) # mean power
  PM[f,11] <- mean(p[1:ppi]) / bodymass[[f]] # relative mean power
  
  PM[f,12] <- (data[[f]]$Total[sj[f] + freq[[f]]*0.05] - data[[f]]$Total[sj[f]])/0.05 # RFD 50 ms
  PM[f,13] <- (data[[f]]$Total[sj[f] + freq[[f]]*0.10] - data[[f]]$Total[sj[f]])/0.10
  PM[f,14] <- (data[[f]]$Total[sj[f] + freq[[f]]*0.15] - data[[f]]$Total[sj[f]])/0.15
  PM[f,15] <- (data[[f]]$Total[sj[f] + freq[[f]]*0.20] - data[[f]]$Total[sj[f]])/0.20
  PM[f,16] <- (data[[f]]$Total[pfi] - data[[f]]$Total[mfi]) / (data[[f]]$Time[pfi] - data[[f]]$Time[mfi]) # RFD from min to max force

  PM[f,17] <- j[(freq[[f]]*0.05)] # impulse at 50, 100, 150, 200, 250 ms after start of movement
  PM[f,18] <- j[(freq[[f]]*0.10)]
  PM[f,19] <- j[(freq[[f]]*0.15)]
  PM[f,20] <- j[(freq[[f]]*0.20)]
  PM[f,21] <- j[length(j)] # impulse at take-off
  
  PM[f,22] <- FT[f] / (data[[f]]$Time[flight[f,1]] - data[[f]]$Time[sj[f]]) # flight time : contraction time
  PM[f,23] <- data[[f]]$Time[flight[f,1]] - data[[f]]$Time[sj[f]] # contraction time
}
rm(f,fz,j,a,v,p,ppi,pfi,mfi)


# Export results -------------------------------------------------------------------------
ex <- readline("Do you want to export the results? (Y/N) ")
if (grepl(ex,'y',ignore.case = T)) {
  tbl <- cbind(JH,PM,ptq)
  rownames(tbl) <- names(data)
  # file name
  fn <- readline("Enter the name of the file you want to save the results to: ")
  # file path
  fp <- dirname(flist[1])

  # write.xlsx(tbl, file = paste0(fp,'/',fn,'.xlsx'))
  write.csv(tbl, file = paste0(fp,'/',fn,'.csv'))
  rm(tbl)
}
