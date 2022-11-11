#' Marginal Means Epi Curve
#' 
#' computeEpiCurve() computes the expected cumulative number of infections/susceptibles per day.
#' 
#' The effects of nonpharmaceutical interventions on an epidemic can be evaluated using this function.  
#' It provides the ability to set a schedule of social distancing/quarantining probabilities, of
#' adaptively reopening society/relaxing social distancing based on a daily hard threshold or by 
#' proportion of the peak daily number of new infections.  Reopening can occur to the entire population 
#' or to a subpopulation based on having low number of contacts.
#' 
#' @param A_mat (sparse) adjacency matrix
#' @param I0 integer.  Number of initial infections in the population
#' @param nDays integer. Number of days to forecast.
#' @param nDaysUntilRecovered integer. Number of days an infective can infect a neighboring susceptible
#' @param probTransPerContact numeric (scalar or vector of length nDays). Daily probability 
#'  that a S-I contact will lead to a new transmission event
#' @param importProb numeric. Daily probability that a susceptible will import the disease from outside 
#'  of the population
#' @param propPeakToReopen numeric. (Adaptive reopening) Quarantining will be relaxed if the number of 
#'  daily new infections reaches propPeakToReopen proportion of the peak.
#' @param minNumberToReopen integer. (Adaptive reopening) Quarantining will be relaxed if the number of 
#'  daily new infections is below minNumberToReopen.  If not NULL, minNumberToReopen will override 
#'  propPeakToReopen.
#' @param minDayToReopen integer. (Adaptive reopening) Reopening will not occur before minDayToReopen. 
#'  This is important if minNumberToReopen is not NULL, sine the epidemic will be below this threshold in 
#'  its early days.
#' @param  gradualReopening integer.  (Adaptive reopening) If not NULL, a Gaussian kernel will be used 
#'  to gradually reopen society such that the quarantining probability will reach 90\% of its current 
#'  level after gradualReopening days.
#' @param nDaysIncreasing integer.  (Adaptive reopening) Number of days of increasing infections 
#'  before reimplementing social distancing/quarantining to secondQuarProb.
#' @param secondQuarProb numeric. (Adaptive reopening) Probability that an individual will be 
#'  quarantined if social distancing needs to be reimplemented.
#' @param quarantineProb numeric vector. Vector of length nDays giving the daily probability that an 
#'  individual will be quarantined.
#' @param quarantineProb_highDegree numeric vector. Vector of length nDays giving the daily 
#'  probability that an individual above highDeg_cutoff will be quarantined.
#' @param highDeg_cutoff integer. The degree cutoff such that anyone with degree higher 
#'  than highDeg_cutoff will follow quarantineProb_highDegree rather than quarantineProb, and if 
#'  adaptive reopening occurs, these high degree individuals will remain at the prespecified 
#'  quarantineProb_highDegree.
#' @param plotProgress logical. 
#' @return A nDays x 5 matrix giving the cumulative number of susceptibles, cumulative number of infected,
#' number of daily infectives (not new infections), and a binary vector showing if social distancing 
#' measures have been relaxed.
#' @export
compute_SEIR <- function(A_mat,
                         I0 = 1,
                         nDays = 30,
                         nDaysUntilRecovered = 11,
                         nDaysSusceptible = 3,
                         probTransPerContact = 0.01,
                         importProb = 1.050940e-06,
                         propPeakToReopen = 0, 
                         minNumberToReopen = NULL,
                         minDayToReopen = Inf, 
                         gradualReopening = 7, 
                         nDaysIncreasing = 7,
                         secondQuarProb = 0.3,
                         quarantineProb = rep(c(0,0.3,0.5),c(6,17,nDays- 17 - 6)),
                         quarantineProb_highDegree = quarantineProb,
                         highDeg_cutoff = 50,
                         plotProgress = TRUE){
  if(class(A_mat) == "igraph") A_mat = as_adjacency_matrix(A)
  
  # Population size
  N = nrow(A_mat)
  
  # Allow for degree-based quarantining
  Degree = rowSums(A_mat)
  highDeg = (Degree > highDeg_cutoff) + 1
  qt = cbind(quarantineProb,quarantineProb_highDegree)
  
  # Too lazy to write out descriptive names
  DI = nDaysUntilRecovered
  DE = nDaysSusceptible
  
  
  # Allow for varying probability of S-I contacts leading to a new infection due to, e.g., face shields
  if(length(probTransPerContact) == 1){
    pTrans = rep(probTransPerContact,nDays)
  }else{
    pTrans = probTransPerContact
  }
  
  # Objects to tally number of infectious individuals at a given time point
  nInf_curr = nInf_cum = integer(nDays)
  
  # Key probability vectors/matrices
  probInf = numeric(N)
  
  # Helpers
  f_s = matrix(1.0,N,DE+DI)
  colSeq = rep(1:(DE+DI),ceiling(nDays/(DE+DI)))
  alphaMat = matrix(0.0,N,DE+DI)
  colSeq_a = rep(1:(DE+DI),ceiling((nDays+1)/(DE+DI)))
  
  # Initial probability anyone in network is infected
  p0 = I0/N
  
  # Initialize counts
  nInf_cum[1] = nInf_curr[1] = sum(I0)
  
  # First day
  probInf[] = p0
  
  # Initialize probability vectors
  alphaMat[,1] = 1.0 - p0
  
  
  # Adaptive objects
  reopening_tracker = adapting_tracker = integer(nDays) #0 if not adapting yet, 1 if adapting
  peakNInf = I0
  nInf_to_reopen = 0
  consecDaysInc = 0
  adapting = reopening = FALSE
  if(!is.null(gradualReopening)){
    s2 = -0.5*gradualReopening^2/log(0.1)
  }
  
  #--- Compare D_E with D_I
  if(DE > DI){
    xtA = probInf%*%A_mat
    
    for(tt in 1:DI){
      # Update \hat{x}_{tj}
      # - Don't need to update x_{sj}
      
      # Update alpha and f_{sj}'s
      # - Don't need to update xtA
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,ifelse(tt>1,diff(nInf_cum[tt + -1:0]),nInf_cum[tt] - I0))
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,
                           main="Current expected number of Infectious",sub=tt)
    }
    
    probInf[] = 0.0
    for(tt in (DI+1):DE){
      # Update \hat{x}_{tj}
      # - Don't need to update x_{sj}
      
      # Update alpha and f_{sj}'s
      # - Don't need to update xtA
      f_s[,colSeq[tt]] = (1-importProb)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = 0
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,diff(nInf_cum[tt + -1:0]))
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,
                           main="Current expected number of Infectious",sub=tt)
    }
    
    for(tt in (DE+1):nDays){
      # Check for adaptivity
      if(reopening){
        consecDaysInc = 
          ifelse(nInf_curr[tt-1] > nInf_curr[tt-2],
                 consecDaysInc + 1,
                 0)
        if(consecDaysInc >= nDaysIncreasing){
          qt[tt:nDays,1] = secondQuarProb
          reopening_tracker[tt:nDays] = 0
          if(!adapting){adapting_tracker[tt:nDays] = 1;adapting = TRUE}
        }
      }
      
      # Update \hat{x}_{tj}
      if(tt == DE+1){
        rs = log(f_s[,colSeq[max(1,tt-DE-DI+1):(tt-DE)]])
      }else{
        rs = rowSums(log(f_s[,colSeq[max(1,tt-DE-DI+1):(tt-DE)]]))
      }
      probInf = alphaMat[,colSeq_a[max(1,tt+1-DE-DI)]]*( 1 - exp(rs) )
      
      # Update alpha
      xtA = probInf%*%A_mat
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,nInf_curr[tt])
      nInf_to_reopen = ifelse(is.null(minNumberToReopen),peakNInf * propPeakToReopen, minNumberToReopen)
      nInf_to_reopen = ifelse(tt >= minDayToReopen, nInf_to_reopen, 0)
      if(!reopening){
        if(nInf_curr[tt] < nInf_to_reopen){
          reopening = TRUE
          reopening_tracker[tt:nDays] = 1
          if(is.null(gradualReopening)){
            qt[tt:nDays,1] = 0
          }else{
            qt[tt:nDays,1] = exp(-0.5*c(0:(nDays-tt))^2/s2)*qt[tt-1]
          }
        }
      }
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,xlab="Time since start",ylab="",
                           main="Current expected number of Infectious",sub=tt)
    }
    
    
    
  }#End: DE > DI
  if(DE < DI){
    xtA = probInf%*%A_mat
    
    for(tt in 1:DE){
      # Update \hat{x}_{tj}
      # - Don't need to update x_{sj}
      
      # Update alpha and f_{sj}'s
      # - Don't need to update xtA
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,ifelse(tt>1,diff(nInf_cum[tt + -1:0]),nInf_cum[tt] - I0))
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,
                           main="Current expected number of Infectious",sub=tt)
    }
    
    rs = numeric(N)
    for(tt in (DE+1):DI){
      # Update \hat{x}_{tj}
      rs = rs + log(f_s[,tt-DE])
      probInf = 1 - (1-p0)*exp(rs)
      
      # Update alpha and f_{sj}'s
      xtA = probInf%*%A_mat
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,diff(nInf_cum[tt + -1:0]))
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,
                           main="Current expected number of Infectious",sub=tt)
    }
    
    for(tt in (DI+1):nDays){
      # Check for adaptivity
      if(reopening){
        consecDaysInc = 
          ifelse(nInf_curr[tt-1] > nInf_curr[tt-2],
                 consecDaysInc + 1,
                 0)
        if(consecDaysInc >= nDaysIncreasing){
          qt[tt:nDays,1] = secondQuarProb
          reopening_tracker[tt:nDays] = 0
          if(!adapting){adapting_tracker[tt:nDays] = 1;adapting = TRUE}
        }
      }
      
      # Update \hat{x}_{tj}
      rs = rowSums(log(f_s[,colSeq[max(1,tt-DE-DI+1):(tt-DE)]]))
      probInf = alphaMat[,colSeq_a[max(1,tt+1-DE-DI)]]*( 1 - exp(rs) )
      
      # Update alpha
      xtA = probInf%*%A_mat
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,nInf_curr[tt])
      nInf_to_reopen = ifelse(is.null(minNumberToReopen),peakNInf * propPeakToReopen, minNumberToReopen)
      nInf_to_reopen = ifelse(tt >= minDayToReopen, nInf_to_reopen, 0)
      if(!reopening){
        if(nInf_curr[tt] < nInf_to_reopen){
          reopening = TRUE
          reopening_tracker[tt:nDays] = 1
          if(is.null(gradualReopening)){
            qt[tt:nDays,1] = 0
          }else{
            qt[tt:nDays,1] = exp(-0.5*c(0:(nDays-tt))^2/s2)*qt[tt-1]
          }
        }
      }
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,xlab="Time since start",ylab="",
                           main="Current expected number of Infectious",sub=tt)
    }
    
    
    
  }#End: DE < DI
  if(near(DE,DI)){
    xtA = probInf%*%A_mat
    
    for(tt in 1:DE){
      # Update \hat{x}_{tj}
      # - Don't need to update x_{sj}
      
      # Update alpha and f_{sj}'s
      # - Don't need to update xtA
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,ifelse(tt>1,diff(nInf_cum[tt + -1:0]),nInf_cum[tt] - I0))
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,
                           main="Current expected number of Infectious",sub=tt)
    }
    
    for(tt in (DE+1):nDays){
      # Check for adaptivity
      if(reopening){
        consecDaysInc = 
          ifelse(nInf_curr[tt-1] > nInf_curr[tt-2],
                 consecDaysInc + 1,
                 0)
        if(consecDaysInc >= nDaysIncreasing){
          qt[tt:nDays,1] = secondQuarProb
          reopening_tracker[tt:nDays] = 0
          if(!adapting){adapting_tracker[tt:nDays] = 1;adapting = TRUE}
        }
      }
      
      # Update \hat{x}_{tj}
      if(tt == DI+1){
        rs = log(f_s[,colSeq[max(1,tt-DE-DI+1):(tt-DE)]])
      }else{
        rs = rowSums(log(f_s[,colSeq[max(1,tt-DE-DI+1):(tt-DE)]]))
      }
      probInf = alphaMat[,colSeq_a[max(1,tt+1-DE-DI)]]*( 1 - exp(rs) )
      
      # Update alpha
      xtA = probInf%*%A_mat
      f_s[,colSeq[tt]] = (1-importProb)*(qt[tt,][highDeg] + (1-qt[tt,][highDeg])*(1-pTrans[tt])^xtA@x)
      alphaMat[,colSeq_a[tt+1]] = alphaMat[,colSeq_a[tt]]*f_s[,colSeq[tt]]
      
      # Update counts
      nInf_curr[tt] = sum(probInf)
      nInf_cum[tt] = N - sum(alphaMat[,colSeq_a[tt+1]])
      
      # Set adaptive values
      peakNInf = max(peakNInf,nInf_curr[tt])
      nInf_to_reopen = ifelse(is.null(minNumberToReopen),peakNInf * propPeakToReopen, minNumberToReopen)
      nInf_to_reopen = ifelse(tt >= minDayToReopen, nInf_to_reopen, 0)
      if(!reopening){
        if(nInf_curr[tt] < nInf_to_reopen){
          reopening = TRUE
          reopening_tracker[tt:nDays] = 1
          if(is.null(gradualReopening)){
            qt[tt:nDays,1] = 0
          }else{
            qt[tt:nDays,1] = exp(-0.5*c(0:(nDays-tt))^2/s2)*qt[tt-1]
          }
        }
      }
      
      if(plotProgress)plot(nInf_curr[1:tt],type='l',lwd=2,xlab="Time since start",ylab="",
                           main="Current expected number of Infectious",sub=tt)
    }
    
    
    
  }#End: DE = DI
  
  return(cbind(nSusc_cum = N - nInf_cum,
               nInf_cum = nInf_cum,
               nInf_daily = nInf_curr,
               reopening = reopening_tracker,
               adapting = adapting_tracker))
}
