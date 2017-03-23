# We show that filtfilt is needed to avoid shifting the signal over. 

signal = sin(1:100)
butter(6, c(13,30)/(128/2), type = 'pass')

signal_fil = signal::filter(bu2, signal)
signal_fil_zero = signal::filtfilt(bu2, signal)

plot(signal, type = "l")  
lines(signal_fil, col = "red", lty = "dotted")
lines(signal_fil_zero, col = "blue",  lty = "dashed")
