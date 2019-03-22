using LinearAlgebra
using WAV
function addnoise(signal,noise,snr)
    SNR(signal,noisy)=20*log10(norm(signal)/norm(signal-noisy))
    S=length(signal)
    N=length(noise)
    @assert N>S "noise should be longer than signal"  
    R = rand(1:1+N-S)
    noise = noise(R:R+S-1)
    noise = noise / norm(noise) * norm(signal) / 10.0^(0.05*snr)
    noisy = signal + noise;
    return noisy,noise
end

signal=wavread("Ideal Binary Mask/matlab/sp10.wav")