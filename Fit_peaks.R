# moorelr/CO2-peaks is licensed under The MIT License
# Copyright 2020 Lowell R. Moore

# NOTE: This script should work as long as the spectra
#   are .txt files in the "../data" directory.

# Note: This script works as a "while" loop which updates a
#   GUI using the "menu" function in base R.
#   - Select "Import spectrum" to load and plot the first
#   spectrum.
#   - Fitting a peak will overlay the data used for the fit
#   in blue and the fit peak curve in red, and append the
#   peak info to an output file. (see "single_peak" function)
#   - Select "Next spectrum" followed by "Import spectrum"
#   to advance to the next spectrum file
#   - Select "Quit" to end the while loop
#   - At the end of the script, the summary data are saved in
#   the working directory as a .csv file

# set working directory ####
dir_save <- getwd()
dir_data <- paste(getwd(), "data", sep = "/")

setwd(dir_data); list.files()
i <- 1      # Start on the first spectrum
n <- sum(grepl(".txt", list.files()))    # how many spectra

print(paste("Processing", n, "spectra in the directory --", getwd()))

# functions used to fit spectra ####

peak <- function(x, center, imax, beta.sq){
  imax*exp(-pi*(x-center)^2/beta.sq)
}

line <- function(x, m, b){(m*x)+b}

single_peak <- function(peak_name, window, wiggling = FALSE, peak_zoom = FALSE){
  
  plot_window <- NULL
  if(peak_zoom){
    plot_window <- window
  }
  
  plot(spectrum$V1, spectrum$V2, cex = 0, main = spectrum_name, xlim = plot_window)
  lines(spectrum$V1, spectrum$V2)
  
  abline(v = window, col = "green", lty = 2)
  
  which.fit <- which(spectrum$V1 >= window[1] & spectrum$V1 <= window[2])                 
  
  xs <- spectrum$V1[which.fit]
  ys <- spectrum$V2[which.fit]
  
  y_start <- mean(ys[1:5])
  y_end <- mean(ys[(length(ys)-5):length(ys)])
  m.i <- (y_end-y_start)/(xs[length(xs)]-xs[1])
  b.i <- y_start-(m.i*xs[1])#as.numeric(quantile(ys, probs = 0.25))
  
  center.a.i <- xs[ys == max(ys)]
  xmid <- (max(xs) + min(xs))/2
  ymid <- (m.i*xmid) + b.i
  imax.a.i <- max(ys) - ymid
  beta.sq.a.i <- 5
  
  if(wiggling == TRUE){
    wiggle <- 0.0001
    imax.a.i <- imax.a.i + (imax.a.i*rnorm(1, sd = wiggle))
    center.a.i <- center.a.i + (center.a.i*rnorm(1, sd = wiggle))
    beta.sq.a.i <- beta.sq.a.i + (beta.sq.a.i*rnorm(1, sd = wiggle))
  }
  
  
  lines(xs, (line(xs, m.i, b.i) + peak(xs, center.a.i, imax.a.i, beta.sq.a.i))
        , col = "blue"
  )
  
  peak.fit <- nls(ys ~ line(xs, m.a, b.a)
                  + peak(xs, center.a, imax.a, beta.sq.a)
                  
                  , start = c(
                    m.a = m.i, b.a = b.i
                    , center.a = center.a.i, imax.a = imax.a.i, beta.sq.a = beta.sq.a.i 
                  )
  )
  
  m.f <- summary(peak.fit)$coefficients[1, 1]
  b.f <- summary(peak.fit)$coefficients[2, 1]
  center.a.f <- summary(peak.fit)$coefficients[3, 1]
  fmax.a.f <- summary(peak.fit)$coefficients[4, 1]
  beta.sq.a.f <- summary(peak.fit)$coefficients[5, 1]
  
  m.f.err <- summary(peak.fit)$coefficients[1, 2]
  b.f.err <- summary(peak.fit)$coefficients[2, 2]
  center.a.f.err <- summary(peak.fit)$coefficients[3, 2]
  fmax.a.f.err <- summary(peak.fit)$coefficients[4, 2]
  beta.sq.a.f.err <- summary(peak.fit)$coefficients[5, 2]
  
  lines(xs, (line(xs, m.f, b.f) + peak(xs, center.a.f, fmax.a.f, beta.sq.a.f))
        , col = "red")
  
  parameters <- cbind.data.frame(spectrum_name, peak_name, m.f, b.f, center.a.f, fmax.a.f, beta.sq.a.f
                                 , m.f.err, b.f.err, center.a.f.err, fmax.a.f.err, beta.sq.a.f.err
  )
  
  return(parameters)
}

# create output dataframe ####

output <- as.data.frame(matrix(nrow = 0, ncol = 12))
colnames(output) <- c("spectrum_name","peak_name", "m.f", "b.f", "center.a.f", "fmax.a.f", "beta.sq.a.f"
                      , "m.f.err", "b.f.err", "center.a.f.err", "fmax.a.f.err", "beta.sq.a.f.err"
)

# Text input mode ####

while(i <= n){
  Sys.sleep(0.5)
  text_input <- menu(c(
      "Fit Ne 1031 peak"
    , "Fit Ne 1042 peak"
    , "Fit Ne 1458 peak"
    , "Fit CO2 1285 peak"
    , "Fit CO2 1388 peak"
    , "Next spectrum"
    , "Import spectrum"
    , "quit"
    , "delete last peak")
    , title = "What to do next?", graphics = TRUE
    )

if(text_input == 6){

  # Go to next spectrum ####

i <- i + 1
print(paste("Now on spectrum", i))

}

if(text_input == 7){
  # Import spectrum ####

spectrum_name <- list.files()[i]

spectrum <- read.delim(spectrum_name, header=FALSE)
plot(spectrum$V1, spectrum$V2, cex = 0, main = spectrum_name)
lines(spectrum$V1, spectrum$V2)

peak_fit <- NA
print(paste("Displaying spectrum", i, "of", n))
}

if(text_input == 1){
  # Fit Ne 1031 peak ####

  for(try in 1:100){
    print(paste("try", as.character(try), sep = " = "))
    try(peak_fit <- single_peak("Ne_1031", c(1020, 1040), wiggling = TRUE)
        , silent = TRUE
    )
    if(is.data.frame(peak_fit)){
      output <- rbind.data.frame(output, peak_fit)
      print(peak_fit)
      peak_fit <- NA
      break
    }
  }  
  
}

if(text_input == 2){
  # Fit Ne 1082 peak ####

  for(try in 1:100){
    print(paste("try", as.character(try), sep = " = "))
    try(peak_fit <- single_peak("Ne_1082", c(1060, 1100), wiggling = TRUE)
        , silent = TRUE
    )
    if(is.data.frame(peak_fit)){
      output <- rbind.data.frame(output, peak_fit)
      print(peak_fit)
      peak_fit <- NA
      break
    }
  } 

}

if(text_input == 3){
  # Fit Ne 1458 peak ####

  for(try in 1:100){
    print(paste("try", as.character(try), sep = " = "))
    try(peak_fit <- single_peak("Ne_1458", c(1448, 1468), wiggling = TRUE)
        , silent = TRUE
    )
    if(is.data.frame(peak_fit)){
      output <- rbind.data.frame(output, peak_fit)
      print(peak_fit)
      peak_fit <- NA
      break
    }
  } 

}

if(text_input == 4){
  # Fit CO2 1285 peak ####

for(try in 1:100){
  print(paste("try", as.character(try), sep = " = "))
  try(peak_fit <- single_peak("CO2_1285", c(1270, 1300), wiggling = TRUE)
      , silent = TRUE
      )
  if(is.data.frame(peak_fit)){
    output <- rbind.data.frame(output, peak_fit)
    print(peak_fit)
    peak_fit <- NA
    break
  }
}

}

if(text_input == 5){
  # Fit CO2 1388 peak ####

for(try in 1:100){
  print(paste("try", as.character(try), sep = " = "))
  try(peak_fit <- single_peak("CO2_1388", c(1373, 1402), wiggling = TRUE)
      , silent = TRUE
      )
  if(is.data.frame(peak_fit)){
    output <- rbind.data.frame(output, peak_fit)
    print(peak_fit)
    peak_fit <- NA
    break
  }
}


}

if(text_input == 8){
  break
}
  
if(text_input == 9){
  # Delete the last output row ####
  output <- output[1:(nrow(output)-1),]
}
# end of while loop ####
}

# export outputs ####

#list.files()
time_stamp <- gsub(Sys.time(), pattern = ":", replacement = "_")
write.csv(output, file = paste(dir_save, "/fit peaks ", time_stamp, ".csv", sep = ""), row.names = FALSE)

# Example spectra ####
if(FALSE){
  pdf(file = paste(dir_save, "Example xl 2 mi 2.pdf", sep = "/")
      , width = 9, height = 4)
  
  par(mfrow = c(1, 3))
  setwd(dir_data); list.files()
  i <- 8      # Start on the first spectrum
  n <- sum(grepl(".txt", list.files()))    # how many spectra
  spectrum_name <- list.files()[i]
  spectrum <- read.delim(spectrum_name, header=FALSE)
  peak_fit <- NA
  print(paste("Displaying spectrum", i, "of", n))
  
  # Panel 1: the whole spectrum
  plot(spectrum$V1, spectrum$V2, type = "l")
  
  # Panel 2: zoom on 1285 peak as points
  plot(spectrum$V1, spectrum$V2, xlim = c(1270, 1300))
  
  # Panel 3: fit of 1285 peak
  for(try in 1:100){
    print(paste("try", as.character(try), sep = " = "))
    try(peak_fit <- single_peak("CO2_1285", c(1270, 1300), wiggling = TRUE, peak_zoom = TRUE)
        , silent = TRUE
    )
    if(is.data.frame(peak_fit)){
      output <- rbind.data.frame(output, peak_fit)
      print(peak_fit)
      peak_fit <- NA
      break
    }
  }
  dev.off()
}
if(FALSE){
  pdf(file = paste(dir_save, "Example peak.pdf", sep = "/")
      , width = 5, height = 4)
  xs <- seq(1, 100, length.out = 100)
  m.i <- 0.1
  b.i <- 10
  center.a.i <- 60
  imax.a.i <- 50
  beta.sq.a.i <- 100
  ys <- line(xs, m.i, b.i) + peak(xs, center.a.i, imax.a.i, beta.sq.a.i)
  plot(xs, ys, type = "l"
       , xlab = "x", ylab = "f(x)")
  points(xs, ys, pch = 19, cex = 0.8, col = rgb(0, 0, 0, 0.3))
  dev.off()
}
