#-----------------------------------------------------------------------------------------
#
# Functions for jump analysis Shiny app
# Robert Schuster (ACU SPRINT)
# September 2022
#
#-----------------------------------------------------------------------------------------

# Find peaks -----------------------------------------------------------------------------
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
importTrial <- function(filepath,filename,jumpType) {
  cols <- c("Time","Z.Left","Z.Right","Left","Right")
  fn <- basename(filename)
  # determine length of header
  skip <- read.csv(filepath, nrows = 20, header = F, blank.lines.skip = F, 
                   col.names = paste0("V",seq_len(max(count.fields(filepath, sep = ',')))), fill = T)
  if (any(skip == "Time", na.rm = T)) {
    skip <- which(skip == 'Time')-1
    
    # load data (without header) and keep only Time, Left and Right columns
    df <- read.csv(filepath, skip = skip)
    if (any(names(df) %in% cols, na.rm = T)) {
      df <- df[,which(names(df) %in% cols)] # keep only left, right and time columns
      colnames(df) <- gsub('Z.','',colnames(df)) # remove 'Z." from column names
      
      df$Total <- df$Left + df$Right
      
      # load header
      header <- read.csv(filepath, nrows = skip-1, header = F)
      # find weight for force normalisation
      if (any(header == "Weight", na.rm = T)) {
        bodymass <- as.numeric(na.omit(header[header == 'Weight',2])) * 9.81
      } else {
        bodymass <- median(df$Total[df$Total >= 300])
      }
      # find frequency
      freq <- as.numeric(na.omit(header[header == 'Frequency',2]))
      if (freq < 1000) {
        message(paste(f,": Sampling frequency =",freq,"\nFrequencies > 1000 Hz are recommended"))
      }
      
      data <- list("df" = df, "bodymass" = bodymass, "freq" = freq, "fn" = fn, "jt" = jumpType)
      return(data)
      
    } else {
      # print warning but continue execution with other files
      warning(paste("The file ", fn, " does not contain the required data and was not processed"))
    }
  } else {
    # print warning but continue execution with other files
    warning(paste("The file", fn, "does not contain the required data and was not processed"))
  }
}

# Determine number of repetitions --------------------------------------------------------
nReps <- function(data) {
  # find peaks (at least 2s apart from each other)
  d <- data$freq*2 # frequency * 2s
  pks <- find_peaks(data$df$Total, m = d)
  # determine a cutoff under which all other peaks cannot be considered IMTP trials
  cutoff <- (max(data$df$Total[pks])-data$bodymass)/2.75 + data$bodymass # 36% of BM normalised max force
  # cutoff <- max(data$df$Total[pks])-250 # within 250 N of max value
  # only keep peaks above cutoff and more than 5s apart
  pks <- pks[which(data$df$Total[pks] > cutoff)]
  d <- data$freq*2.5 # frequency * 2.5s
  if (any(diff(pks) < d)) {
    pks <- pks[-(which(diff(pks) < d)+1)]
  }
  
  # If trial contains more than one rep, cut trial into separate reps
  # determine start and end of rep
  for (p in 1:length(pks)) {
    # define arbitrary period to look for start of rep
    if (p > 1) { # if not first rep
      ds <- min(data$freq*3, (pks[p]-pks[p-1])/2) # start of period = halfway between two adjacent peaks or 3S before peak
    } else {
      ds <- min(data$freq*3, pks[p]-1) # else start of period = either start of trial or 3s before peak
    }
    
    # define end of rep
    if (p < length(pks)) { # if not last rep
      de <- min(data$freq*2, (pks[p+1]-pks[p])/2) # end of rep = halfway between two adjacent peaks or 2s after peak
    } else {
      de <- min(data$freq*2, length(data$df$Total)-pks[p]) # else end of rep = either end of trial or 2s after peak
    }
    w = (pks[p]-ds):(pks[p]+de)
    
    if (length(pks) > 1) {
      n <- paste0(data$fn,'_',p)
      data[[n]] <- data$df[w,]
      
      data$fmaxi[p] <- (pks[p] - w[1]) + 1
      
      if (p == length(pks)) {
        data$df <- NULL
      }
    } else {
      data[[data$fn]] <- data$df[w,]
      data$fmaxi <- (pks - w[1]) + 1
      
      data$df <- NULL
    }
  }
  return(data)
}


# Determine start and end of flight ------------------------------------------------------
flight <- function(data, thl) {
  reps <- names(data)[grep(data$fn,names(data),fixed = T)]
  sj <- numeric(length(reps))
  flight <- matrix(0,length(reps),2)
  colnames(flight) <- c('start','end')
  for (r in 1:length(reps)) {
    rn <- reps[r]
    # # threshold = 5 SD of 1s weighing period before flight or smallest SD or 1s rolling window
    # th <- sd(data[[rn]]$Total[1:(data$freq*1)])*5
    # th <- min(movstd(data[[rn]]$Total, (data$freq*1)), na.rm = T)*5
    th <- min(movstd(data[[rn]]$Total, (data$freq*thl)), na.rm = T)*5
    
    # flight thresholds = 5N start, 20N end
    sft <- 5
    eft <- 20
    
    # start and end of flight = first and last instance when force < flight thresholds
    zero <- find_flight(data[[rn]]$Total[1:data$fmaxi[r]], sft)
    flight[r,2] <- min(which(data[[rn]]$Total[zero$start:length(data[[rn]]$Total)] > eft)) + zero$start
    ep <- min(which(data[[rn]]$Total[1:flight[r,2]] == max(data[[rn]]$Total[1:flight[r,2]])))
    flight[r,1] <- min(which(data[[rn]]$Total[ep:flight[r,2]] < sft)) + ep
    em <- min(which(data[[rn]]$Total[(ep-data$freq*1):ep] == min(data[[rn]]$Total[(ep-data$freq*1):ep]))) + (ep-data$freq*1)
    
    # start of jump movement (sj)
    if (data$jt == 'Countermovement jump') {
      if (any(data[[rn]]$Total[1:em] < (data$bodymass - th))) { # there is a countermovement
        sj[r] <- max(which(data[[rn]]$Total[1:em] > (data$bodymass - th)))
        # sj[r] <- max(which(data[[rn]]$Total[1:em] > (data$bodymass - th))) - (data$freq*0.03)
      } else { # remove jump
        sj[r] <- NULL
        flight[r,] <- NULL
        data[[rn]] <- NULL
      }
    } else if (data$jt == 'Squat jump') {
      sj[r] <- max(which(data[[rn]]$Total[1:ep] < (data$bodymass + th)))
      # sj[r] <- max(which(data[[rn]]$Total[1:ep] < (data$bodymass + th))) - (data$freq*0.03)
    }
  }
  data$flight <- flight
  data$sj <- sj
  return(data)
}


# Check for indicators of poor trial -----------------------------------------------------
qualityCheck <- function(data) {
  reps <- names(data)[grep(data$fn,names(data),fixed = T)]
  data$warn <- list()
  for (r in 1:length(reps)) {
    rn <- reps[r]
    sf <- data$flight[r,1]
    th <- min(movstd(data[[rn]]$Total, (data$freq*1)), na.rm = T)*5
    i <- 1
    
    # unstable weighing period before jump (change in force > 50 N)
    if (any(diff(data[[rn]]$Total[1:sf]) > 50)) {
      data$warn[[rn]][[i]] <- "This rep does not have a stable weighing period"
      i <- i + 1
    }
    # countermovement before jump
    if (data$jt == 'Squat jump' && any(data[[rn]]$Total[(data$sj[r] - data$freq*0.25):data$sj[r]] < (data$bodymass - th))) {
      data$warn[[rn]][[i]] <- "This rep contains a countermovement prior to the squat jump"
      i <- i + 1
    }
  }
  return(data)
}


# Extract performance metrics ------------------------------------------------------------
perfMetrics <- function(data) {
  reps <- names(data)[grep(data$fn,names(data),fixed = T)]
  pm <-  matrix(0,length(reps),29)
  colnames(pm) <- c('JH_ft','JH_J','JH_di','v t-o',
                    'peak force','rel peak force','peak velocity','peak power','rel peak power',
                    'mean force','rel mean force','mean velocity','mean power','rel mean power',
                    'RFD-50','RFD-100','RFD-150','RFD-200','RFD-MM',
                    'J-50','J-100','J-150','J-200','J-T',
                    'FT:CT','Contraction T','Concentric T','Eccentric T','CM Displacement')
  rownames(pm) <- reps
  
  for (r in 1:length(reps)) {
    rn <- reps[r]
    sj <- data$sj[r]
    # jump height
    # flight time (FT)
    FT <- data[[rn]]$Time[data$flight[r,2]] - data[[rn]]$Time[data$flight[r,1]]
    u <- 9.81 * (FT/2) # velocity at take-off
    pm[r,1] <- (0^2 - u^2) / (2*-9.81)
    # impulse (J)
    Fr <- data[[rn]]$Total[sj:data$flight[r,1]] - data$bodymass # resultant force between start of movement and take-off
    J <- sum(diff(data[[rn]]$Time[sj:data$flight[r,1]]) * (head(Fr,-1) + tail(Fr,-1)))/2 # impulse (trapz integration of resultant force)
    pm[r,2] <- ((J/(data$bodymass/9.81))^2)/(2*9.81)
    # double integration (di)
    a <- (data[[rn]]$Total[sj:data$flight[r,2]] / (data$bodymass / 9.81)) - 9.81 # acceleration between start of movement and take-off
    v <- cumsum(c(0,(a[1:(length(a)-1)] + a[2:length(a)])/2 * diff(data[[rn]]$Time[sj:data$flight[r,2]]))) # velocity (trapz integration of acceleration)
    vp <- v[v >= 0] # positive velocity (equivalent to concentric only)
    d <- cumsum(c(0,(v[1:(length(v)-1)] + v[2:length(v)])/2 * diff(data[[rn]]$Time[sj:data$flight[r,2]]))) # displacement (trapz integration of velocity)
    pm[r,3] <- max(d) - d[length(d)]
    # stat of concentric (sc) & start of positive acceleration (sp)
    if (data$jt == 'Countermovement jump') {
      ft <- data$flight[r,2] - data$flight[r,1]
      sc <- min(which(d[1:(length(d)-ft)] == min(d[1:(length(d)-ft)]))) + sj # bottom most position of countermovement
      sp <- min(which(a >= 0)) + sj
    } else if (data$jt == 'Squat jump') {
      sc <- sj
      sp <- sj
    }
    
    # impulse
    fz <- data[[rn]]$Total[sc:data$flight[r,1]] # Fz = force between start of concentric and take-off
    j <- cumsum(c(0,(fz[1:(length(fz)-1)] + fz[2:length(fz)])/2 * diff(data[[rn]]$Time[sc:data$flight[r,1]]))) # impulse (trapz integration of force)
    # power
    p <- data[[rn]]$Total[sj:data$flight[r,2]] * v # power = Fz * velocity
    # peak power index
    ppi <- min(which(p == max(p)))
    # peak & min force indices
    pfi <- min(which(data[[rn]]$Total[sj:data$flight[r,1]] == max(data[[rn]]$Total[sj:data$flight[r,1]]))) + sj
    mfi <- min(which(data[[rn]]$Total[sj:pfi] == min(data[[rn]]$Total[sj:pfi]))) + sj
    
    pm[r,4] <- v[data$flight[r,1]-sj] # take-off velocity
    
    pm[r,5] <- max(data[[rn]]$Total[sj:data$flight[r,1]]) # peak force
    pm[r,6] <- max(data[[rn]]$Total[sj:data$flight[r,1]]) / (data$bodymass/9.81) # relative peak force
    pm[r,7] <- max(v) # peak velocity
    pm[r,8] <- max(p) # peak power
    pm[r,9] <- max(p) / (data$bodymass/9.81) # relative peak power
    
    pm[r,10] <- mean(data[[rn]]$Total[sj:data$flight[r,1]]) # mean force
    pm[r,11] <- mean(data[[rn]]$Total[sj:data$flight[r,1]]) / (data$bodymass/9.81) # relative mean force
    pm[r,12] <- mean(vp[1:max(which(vp == max(vp)))]) # mean velocity
    pm[r,13] <- mean(p[1:ppi]) # mean power
    pm[r,14] <- mean(p[1:ppi]) / (data$bodymass/9.81) # relative mean power
    
    pm[r,15] <- (data[[rn]]$Total[sp + data$freq*0.05] - data[[rn]]$Total[sp])/0.05 # RFD 50 ms
    pm[r,16] <- (data[[rn]]$Total[sp + data$freq*0.10] - data[[rn]]$Total[sp])/0.10
    pm[r,17] <- (data[[rn]]$Total[sp + data$freq*0.15] - data[[rn]]$Total[sp])/0.15
    pm[r,18] <- (data[[rn]]$Total[sp + data$freq*0.20] - data[[rn]]$Total[sp])/0.20
    pm[r,19] <- (data[[rn]]$Total[pfi] - data[[rn]]$Total[mfi]) / (data[[rn]]$Time[pfi] - data[[rn]]$Time[mfi]) # RFD from min to max force
    
    pm[r,20] <- j[data$freq*0.05] # impulse at 50, 100, 150, 200 ms after start of concentric
    pm[r,21] <- j[data$freq*0.10]
    pm[r,22] <- j[data$freq*0.15]
    pm[r,23] <- j[data$freq*0.20]
    pm[r,24] <- j[length(j)] # impulse at take-off
    
    pm[r,25] <- FT / (data[[rn]]$Time[data$flight[r,1]] - data[[rn]]$Time[data$sj[r]]) # FT:CT
    pm[r,26] <- data[[rn]]$Time[data$flight[r,1]] - data[[rn]]$Time[data$sj[r]] # contraction time
    if (data$jt == 'Countermovement jump') {
      pm[r,27] <- data[[rn]]$Time[sc] - data[[rn]]$Time[data$sj[r]] # eccentric time
      pm[r,28] <- data[[rn]]$Time[data$flight[r,1]] - data[[rn]]$Time[sc] # concentric time
      pm[r,29] <- min(d[1:(length(d)-(data$flight[r,2] - data$flight[r,1]))]) # bottom most position of countermovement
    } else if (data$jt == 'Squat jump') {
      pm[r,27] <- NA # eccentric time
      pm[r,28] <- NA # concentric time
      pm[r,29] <- NA # bottom most position of countermovement
    }
  }
  data$pm <- round(pm,2)
  return(data)
}

