using LinearAlgebra
using WAV
function addnoise(signal,noise,snr)
    SNR(signal,noisy)=20*log10(norm(signal)/norm(signal-noisy))
    S=length(signal)
    N=length(noise)
    @assert N>S "noise should be longer than signal"  
    R = rand(1:1+N-S)
    noise = noise[R:R+S-1]
    noise = noise / norm(noise) * norm(signal) / 10.0^(0.05*snr)
    noisy = signal + noise;
    return noisy,noise
end

function vec2frames(data::Array,samples_frame::Int,hops::Int)
    frames=[]
    while length(data) >= samples_frame
        frame=collect(Iterators.take(data,samples_frame))
        data=Iterators.drop(data,hops)
        push!(frames,frame)
    end
    return hcat(frames...)
end

signal=vec(first(wavread("Ideal Binary Mask/matlab/sp10.wav")))
noise=vec(first(wavread("Ideal Binary Mask/matlab/ssn.wav")))
mixed,noise=addnoise(signal,noise,-5)