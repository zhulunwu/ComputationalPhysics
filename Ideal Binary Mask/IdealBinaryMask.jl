function addnoise(signal,noise,snr)
    SNR = @(signal,noisy)( 20*log10(norm(signal)/norm(signal-noisy)) );    
end