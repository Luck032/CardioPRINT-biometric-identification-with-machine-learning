if (1) {
  
  rm(list=ls())

  setwd("cardiprint_signals")
  path = "cardiprint_signals"
  fs = 2000 # acquisition sampling frequence

  data <- read.csv("timestamps_with_neutral.csv")
  Examinee <- data$ID
  for (examinee in Examinee)
  {
    sequence <- 1 
    
    sig <- read.table(paste(examinee,".txt", sep = ''), header = FALSE)

    T1 <- data$T1[which.max(data$ID==examinee)]*fs # start of resting phase
    T2 <- data$T2[which.max(data$ID==examinee)]*fs # end of resting phase
    T3 <- data$T3[which.max(data$ID==examinee)]*fs # start of stress phase
    T4 <- data$T4[which.max(data$ID==examinee)]*fs # end of stress phase
    T5 <- data$T5[which.max(data$ID==examinee)]*fs # start of relaxation phase
    T6 <- data$T6[which.max(data$ID==examinee)]*fs # end of relaxation phase
    
    colnames(sig)[c(1,2)] <- c("ECG","ICG")
    sig <- sig[,c(1,2)]
    sig$time <- (1:length(sig$ECG))/fs
    
    library(signal)
    
    filt1 <- butter(4,c(1,40)/(fs/2), type = "pass")
    filt2 <- butter(4,c(0.5,40)/(fs/2), type = "pass")
    
    sig$fECG <- filtfilt(filt1, sig$ECG) # highpass 1Hz, lowpass 40Hz
    sig$fICG <- filtfilt(filt2, sig$ICG) # highpass  0.5Hz, lowpass 40Hz
    
    # Pan-Tompkins method of QRS peak detection
    
    i <- 1
    L <- length(sig$time)
    N <- 150 * fs/1000 # movingwaveform window of 150 ms
    sig$fECGmov <- vector("integer", L)
    sig$fECGmov[1:(L - 1)] <- sig$fECG[2:L]
    sig$ECGdiff <- (sig$fECGmov - sig$fECG)/fs
    sig$quad <- sig$ECGdiff^2
    sig$MovWav <- vector("integer", L)
    
    while (i < (N + 1))
    {
      sig$MovWav[N:L] <- sig$MovWav[N:L] + sig$quad[i:(L - N + i)]
      i <- i + 1
    }

    j <- 0
    maximum_corrected <- 0
    while (j<10)
    {
      maximum_corrected <- maximum_corrected + max(sig$MovWav[round(j*L/10):(round(j*L/10)+5000)])
      j <- j + 1
    }
    
    maximum_corrected <- maximum_corrected/10
    
    threshold <- 0.3 * maximum_corrected # setting threshold based on maximal value in entire signal
    blanking <- (150 + 50) * fs/1000   # 50 ms longer than mowing waveform window to ensure there is no R wave double-counting
    i <- blanking                      # offset to skip R wave in the begining of the signal, since P-QRS-T segment would not be complete and therefore not needed for this analysis
    ii <- 1
    Rpeak <- 0
    
    while (i < L)                      # while loop makes sure that there iterations dont pass through maximum samples of data
    {
      if (sig$MovWav[i] > threshold)   # if current signal value is higher than set threshold it means that QRS peak occured
      { 
        Rpeak[ii] <- i                  # in variable Rpeak we save location of peaks
        i <- i + blanking
        ii <- ii + 1
      }
      
      else
      {
        i <- i + 1
      }
    }

    i <- 1
    M <- length(Rpeak)
    Cpeak <- 0
    
    while (i < M)
    {
      maxim <- which.max(sig$fECG[(Rpeak[i] - 0.03*fs):(Rpeak[i] + 0.03*fs)])
      Rpeak[i] <- Rpeak[i] + maxim - 1 - 0.03*fs
      Cpeak[i] <- Rpeak[i] + which.max(sig$fICG[Rpeak[i]:Rpeak[i+1]]) - 1       # location of C peak calculated as maximum value between two R peaks
      RC <- 1000*(Cpeak[i] - Rpeak[i])/fs
      i <- i + 1
    }
    
    print(mean(Rpeak[2:length(Rpeak)]/fs-Rpeak[1:(length(Rpeak)-1)]/fs))
    print(sd(Rpeak[2:length(Rpeak)]/fs-Rpeak[1:(length(Rpeak)-1)]/fs))

    while (sequence < 4) {
      
      if (sequence == 1)
      {
        dRpeak <- Rpeak[Rpeak > T1]
        dRpeak <- dRpeak[dRpeak < T2]
        dCpeak <- Cpeak[Rpeak > T1]
        dCpeak <- dCpeak[1:length(dRpeak)]
      }
      
      else
      {
        if (sequence == 2)
        {
          dRpeak <- Rpeak[Rpeak > T3]
          dRpeak <- dRpeak[dRpeak < T4]
          dCpeak <- Cpeak[Rpeak > T3]
          dCpeak <- dCpeak[1:length(dRpeak)]
        }
        
        else
        {
          # if (sequence == 3)
          # {
            dRpeak <- Rpeak[Rpeak > T5]
            dRpeak <- dRpeak[dRpeak < T6]
            dCpeak <- Cpeak[Rpeak > T5]
            dCpeak <- dCpeak[1:length(dRpeak)]
          # }
        }
      }
      
      RCmean <- round(mean(dCpeak - dRpeak))
      
      i <- 1
      segment <- 0.75*fs                                      # length of each segment which will be extracted is 750 ms
      N <- length(dRpeak)                                     # number of peaks in one interval
      N1 <- 10                                                # number of PQRST and BCX segments of ECG and ICG signals, which will be used for averaging
      N2 <- c(2,round(N/3) - 4,round(2*N/3) - 4,N-11)         # orientation of 4 sequences of 10 peaks in order to get more features from every single interval
      PQRST <- matrix(0, segment, 4*N1)
      BCX <- matrix(0, segment, 4*N1)
      RR <- matrix(0, 1, 4*(N1 - 1))
      CC <- matrix(0, 1, 4*(N1 - 1))
      RC <- matrix(0, 1, 4*(N1 - 1))
      time <- (1:segment)/fs
      
      ECG1 <- 0                                               # these are variables for 2x4 different segments in interval, where 10 PQRST intervals will be saved
      ECG2 <- 0                                               
      ECG3 <- 0
      ECG4 <- 0
      ICG1 <- 0                                              
      ICG2 <- 0                                               
      ICG3 <- 0
      ICG4 <- 0
      
      while (i < (N1 + 1))
      {
        PQRST[,i] <- sig$fECG[(dRpeak[N2[1] + i - 1] - 0.25*fs + 1):(dRpeak[N2[1] + i - 1] + 0.5*fs)]
        BCX[,i] <- sig$fICG[(dCpeak[N2[1] + i - 1] - 0.25*fs + 1):(dCpeak[N2[1] + i - 1] + 0.5*fs)]
        RR[i] <- (dRpeak[N2[1] + i] - dRpeak[N2[1] + i - 1])/fs
        CC[i] <- (dCpeak[N2[1] + i] - dCpeak[N2[1] + i - 1])/fs
        RC[i] <- (dCpeak[N2[1] + i] - dRpeak[N2[1] + i])/fs
        ECG1 <- ECG1 + PQRST[,i]
        ICG1 <- ICG1 + BCX[,i]
        
        PQRST[,10 + i] <- sig$fECG[(dRpeak[N2[2] + i - 1] - 0.25*fs + 1):(dRpeak[N2[2] + i - 1] + 0.5*fs)]
        BCX[,10 + i] <- sig$fICG[(dCpeak[N2[2] + i - 1] - 0.25*fs + 1):(dCpeak[N2[2] + i - 1] + 0.5*fs)]
        RR[10 + i] <- (dRpeak[N2[2] + i] - dRpeak[N2[2] + i - 1])/fs
        CC[10 + i] <- (dCpeak[N2[2] + i] - dCpeak[N2[2] + i - 1])/fs
        RC[10 + i] <- (dCpeak[N2[2] + i] - dRpeak[N2[2] + i])/fs
        ECG2 <- ECG2 + PQRST[,10 + i]
        ICG2 <- ICG2 + BCX[,10 + i]
        
        PQRST[,20 + i] <- sig$fECG[(dRpeak[N2[3] + i - 1] - 0.25*fs + 1):(dRpeak[N2[3] + i - 1] + 0.5*fs)]
        BCX[,20 + i] <- sig$fICG[(dCpeak[N2[3] + i - 1] - 0.25*fs + 1):(dCpeak[N2[3] + i - 1] + 0.5*fs)]
        RR[20 + i] <- (dRpeak[N2[3] + i] - dRpeak[N2[3] + i - 1])/fs
        CC[20 + i] <- (dCpeak[N2[3] + i] - dCpeak[N2[3] + i - 1])/fs
        RC[20 + i] <- (dCpeak[N2[3] + i] - dRpeak[N2[3] + i])/fs
        ECG3 <- ECG3 + PQRST[,20 + i]
        ICG3 <- ICG3 + BCX[,20 + i]
        
        PQRST[,30 + i] <- sig$fECG[(dRpeak[N2[4] + i - 1] - 0.25*fs + 1):(dRpeak[N2[4] + i - 1] + 0.5*fs)]
        BCX[,30 + i] <- sig$fICG[(dCpeak[N2[4] + i - 1] - 0.25*fs + 1):(dCpeak[N2[4] + i - 1] + 0.5*fs)]
        RR[30 + i] <- (dRpeak[N2[4] + i] - dRpeak[N2[4] + i - 1])/fs
        CC[30 + i] <- (dCpeak[N2[4] + i] - dCpeak[N2[4] + i - 1])/fs
        RC[30 + i] <- (dCpeak[N2[4] + i] - dRpeak[N2[4] + i])/fs
        ECG4 <- ECG4 + PQRST[,30 + i]
        ICG4 <- ICG4 + BCX[,30 + i]
        
        i <- i + 1
      }
      
      par(new = FALSE)
      
      ECG1 <- ECG1/N1
      ECG2 <- ECG2/N1
      ECG3 <- ECG3/N1
      ECG4 <- ECG4/N1
      ICG1 <- ICG1/N1
      ICG2 <- ICG2/N1
      ICG3 <- ICG3/N1
      ICG4 <- ICG4/N1
      RR1 <- round(sum(RR[1:10])/10,4)
      RR2 <- round(sum(RR[11:20])/10,4)
      RR3 <- round(sum(RR[21:30])/10,4)
      RR4 <- round(sum(RR[31:40])/10,4)
      CC1 <- round(sum(CC[1:10])/10,4)
      CC2 <- round(sum(CC[11:20])/10,4)
      CC3 <- round(sum(CC[21:30])/10,4)
      CC4 <- round(sum(CC[31:40])/10,4)
      RC1 <- round(sum(RC[1:10])/10,4)
      RC2 <- round(sum(RC[11:20])/10,4)
      RC3 <- round(sum(RC[21:30])/10,4)
      RC4 <- round(sum(RC[31:40])/10,4)
      rm(RR,CC,RC)
  
      # In next segment we will define variables which are used to define waves in our signals
      ii <- 1
      R <- 0.25 # in seconds
      Q <- vector(mode = 'integer',length = 4)
      S <- vector(mode = 'integer',length = 4)
      T <- vector(mode = 'integer',length = 4)
      Tstart <- vector(mode = 'integer',length = 4)
      Tend <- vector(mode = 'integer',length = 4)
      Tleft90 <- vector(mode = 'integer',length = 4)
      Tright90 <- vector(mode = 'integer',length = 4)
      ECGcrest <- vector(mode = 'integer',length = 4)
      ECGQRScrest <- vector(mode = 'integer',length = 4)
      ECGTcrest <- vector(mode = 'integer',length = 4)
      C <- 0.25 # in seconds
      B <- vector(mode = 'integer',length = 4)
      X <- vector(mode = 'integer',length = 4)
      ICGcrest <- vector(mode = 'integer',length = 4)

      
      library(dplyr)
      
      # we will use while loop to extract features using all 4 segments of 10 intervals
      while (ii < 5) {
        
        ECG <- case_when(
          ii == 1 ~ ECG1,
          ii == 2 ~ ECG2,
          ii == 3 ~ ECG3,
          ii == 4 ~ ECG4
        )
        
        ICG <- case_when(
          ii == 1 ~ ICG1,
          ii == 2 ~ ICG2,
          ii == 3 ~ ICG3,
          ii == 4 ~ ICG4
        )
        
        RR <- case_when(
          ii == 1 ~ RR1,
          ii == 2 ~ RR2,
          ii == 3 ~ RR3,
          ii == 4 ~ RR4
        )
        
        CC <- case_when(
          ii == 1 ~ CC1,
          ii == 2 ~ CC2,
          ii == 3 ~ CC3,
          ii == 4 ~ CC4
        )
        
        RC <- case_when(
          ii == 1 ~ RC1,
          ii == 2 ~ RC2,
          ii == 3 ~ RC3,
          ii == 4 ~ RC4
        )
        
        dECG <- (ECG[2:(fs*0.75)] - ECG[1:(fs*0.75-1)])/fs          # first derivative of ECG
        
        i <- 0.01*fs                                                # these i & j incremental variables will be used to determine local minimums right and left from R peak, which will represent Q and S points
        j <- 0.01*fs                                                # 10ms offset will be used, as we don't expect Q and S points to be within that period from R peak
        
        while (dECG[0.25*fs - i] > 0)
        {
          i <- i + 1
        }
        
        while(dECG[0.25*fs + j] < 0)
        {
          j <- j + 1
        }
        
        Q[ii] <- 0.25*fs - i
        S[ii] <- 0.25*fs + j
      
        
        T[ii] <- which.max(ECG[(S[ii] + fs*50/1000):(segment - fs*200/1000)])
        T[ii] <- T[ii] + S[ii] + fs/20 - 1

        ECGmov2 <- vector("integer", segment-2)
        d2ECG <- vector("integer", segment-2)
        ECGmov2[1:(segment - 2)] <- dECG[2:(segment-1)]
        d2ECG <- ECGmov2 - dECG[1:(segment-2)]
        
        # Aright
        i <- 1
        while (d2ECG[T[ii] + i + 10] < 0)
        {
          i <- i + 1
        }
        
        Aright <- T[ii] + i + 10
        Bright <- Aright + 200
        
        # Aleft
        i <- 1
        while (d2ECG[T[ii] - i - 10] < 0)
        {
          i <- i + 1
        }
        
        Aleft <- T[ii] - i - 10
        Bleft <- S[ii] + 100
        
        #desni trapez
        i <- 1
        Pright <- vector("integer", Bright-Aright-1)
        while (Aright + i < Bright)
        {
          Pright[i] <- 0.5*(ECG[Aright]-ECG[Aright+i])*(2*Bright-2*Aright-i)
          i<-i+1
        }
        
        Tend[ii] <- Aright + which.max(Pright)
        
        
        
        #levi trapez
        i <- 1
        Pleft <- vector("integer", Aleft-Bleft-1)
        while (Aleft - i > Bleft)
        {
          Pleft[i] <- abs(0.5*(ECG[Aleft]-ECG[Aleft-i])*(2*Bleft-2*Aleft+i))
          i<-i+1
        } 
      
        Tstart[ii] <- Aleft - which.max(Pleft)
        
        i <- 1
        
        while (ECG[Tend[ii] - i] < (ECG[Tend[ii]] + 0.9*(ECG[T[ii]] - ECG[Tend[ii]])))
        {
          i <- i + 1
        }
        
        Tright90[ii] <- Tend[ii] - i + 1                                    #  Right sided point on T wave where it's value is 90% of T wave relative amplitude
        
        i <- 1
        
        while (ECG[Tstart[ii] + i] < (ECG[Tstart[ii]] + 0.9*(ECG[T[ii]] - ECG[Tstart[ii]])))
        {
          i <- i + 1
        }
        
        Tleft90[ii] <- Tstart[ii] + i - 1 
        
        # Crest faktor QRS talasa
        
        i <- Q[ii]
        j <- 0
        CREST <- 0
        
        while (i < (S[ii] + 1))
        {
          amplitude <- (ECG[i] - ECG[Q[ii]])^2
          CREST <- CREST + amplitude
          j <- j+1
          i <- i+1
        }
        
        ECGQRScrest[ii] <- (ECG[T[ii]] - ECG[Q[ii]])/sqrt(CREST/j)
        
        
        
        # Crest faktor T wave
        
        i <- Tstart[ii]
        j <- 0
        CREST <- 0
        
        while (i < (Tend[ii] + 1))
        {
          amplitude <- (ECG[i] - ECG[Tstart[ii]])^2
          CREST <- CREST + amplitude
          j <- j+1
          i <- i+1
        }
        
        ECGTcrest[ii] <- (ECG[T[ii]] - ECG[Tstart[ii]])/sqrt(CREST/j)
        
        # analysis of BCX segments of ICG
        
        RC <- RC*1000 # convert to miliseconds from seconds
        
        # in order to find X point, we will need first derivative of ICG segment to find local minimums
        
        ICGmov <- vector("integer", segment-1)
        dICG <- vector("integer", segment-1)
        dICGmov <- vector("integer", segment-2)
        d2ICG <- vector("integer", segment-2)
        ICGmov[1:(segment - 1)] <- ICG[2:segment]
        dICG <- ICGmov - ICG[1:(segment - 1)]
        dICGmov[1:(segment - 2)] <- dICG[2:(segment-1)]
        d2ICG <- dICGmov - dICG[1:(segment - 2)]
        B[ii] <- which.max(d2ICG[((C-0.080)*fs):(C*fs)]) + (C-0.080)*fs
        
        potentialX <- vector("integer",1)
        
        i <- 1
        j <- 1
        while ((B[ii] + 230*fs/1000 + i) < (B[ii] + 400*fs/1000))
        {
          if (dICG[(B[ii] + 230*fs/1000 + i)] < 0 && dICG[(B[ii] + 230*fs/1000 + i + 1)] > 0) 
          { 
            potentialX[j] <- B[ii] + 230*fs/1000 + i
            j <- j + 1
          } 
          
          i <- i + 1 
        }
        
        if (potentialX[1] == 0) {
          i <- 1
          j <- 1
          while ((B[ii] + 200*fs/1000 + i) < (B[ii] + 400*fs/1000))
          {
            if (dICG[(B[ii] + 200*fs/1000 + i)] < 0 && dICG[(B[ii] + 200*fs/1000 + i + 1)] > 0) 
            { 
              potentialX[j] <- B[ii] + 200*fs/1000 + i
              j <- j + 1
            } 
            
            i <- i + 1 
          }
          
          }
        
        X[ii] <- potentialX[which.min(ICG[potentialX])]
        
        # Crest faktor BC talasa, sa pocetkom u B tacki a kraj u parametru oznacenom kao pik1
        
        i <- B[ii]
        j <- 0
        CREST <- 0
        
        while (i < (X[ii] + 1))
        {
          amp <- (ICG[i] - ICG[B[ii]])^2
          CREST <- CREST + amp
          j <- j+1
          i <- i+1
        }
        
        ICGcrest[ii] <- (ICG[C*fs] - ICG[B[ii]])/sqrt(CREST/j)
        
        # ploting both PQRST and BCX segments
        
        xinterval <- c(1/fs,1)
        yinterval <- c(1.1*min(min(ECG,ICG)),1.1*max(max(ECG,ICG)))

        # calculating amlpitude and temporal features of ECG and ICG
        
        QRS_int <- (S[ii] - Q[ii])*1000/fs # in miliseconds
        RS_amp <- ECG[R*fs] - ECG[S[ii]]
        RQ_amp <- ECG[R*fs] - ECG[Q[ii]]
        RT_amp <- ECG[R*fs] - ECG[T[ii]]
        
        
        T_int <- (Tend[ii] - Tstart[ii])*1000/(fs*sqrt(RR))
        QT_int <- (Tend[ii] - Q[ii])*1000/(fs*sqrt(RR))
        ST_int <- (Tend[ii] - S[ii])*1000/(fs*sqrt(RR))
        TT1_amp <- ECG[T[ii]] - ECG[Tstart[ii]]
        TT2_amp <- ECG[T[ii]] - ECG[Tend[ii]]
        
        RC_int <- RC*1000/fs # left-chamber performance index
        RC_int2 <- RC_int/(sqrt(RR))
        BX_int <- (X[ii] - B[ii])*1000/fs 
        BX_int2 <- BX_int/sqrt(RR)
        CB_amp <- ICG[C*fs] - ICG[B[ii]]
        CX_amp <- ICG[C*fs] - ICG[X[ii]]
        
        features <- round(c(QRS_int, T_int, QT_int, RC_int, RC_int2, BX_int, BX_int2, RS_amp, RQ_amp, TT1_amp, TT2_amp, CB_amp, CX_amp,ECGTcrest[ii], ICGcrest[ii], 1000*RR, 1000*CC),2)
        print(sequence)
        features <- t(c(paste(examinee,'-',sequence), features))
        features <- as.data.frame(features)
        colnames(features) <- c('examinee', 'QRS_int', 'T_int', 'QT_int', 'RC_int', 'RC_int2', 'BX_int'
                                , 'BX_int2', 'RS_amp', 'RQ_amp', 'TT1_amp', 'TT2_amp', 'CB_amp', 'CX_amp', 'ECGTcrest', 'ICGcrest', 'RR', 'CC')
        print(features$examinee)
        if(file.exists('/Fetureset/Features.txt')) {
  
          features2 <- read.table('/Fetureset/Features.txt', sep = '\t', dec = ",")
          rownames(features) <- as.character(nrow(features2)+1)
          features <- rbind(features2,features)
  
        }
  
        write.table(features,'/Fetureset/Features.txt', sep = '\t', dec = ",")
        library(writexl)
        write.csv(x = features, file = "/Fetureset/Features.csv", row.names = FALSE)


        
        ii <- ii + 1
        
      }

      sequence <- sequence + 1
      
    }   # end of while loop for switching between phases
    
    print(examinee)
  }
}


