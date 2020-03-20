using GR
function rope_wave()
    T=0.8
    ρ=0.2
    v=sqrt(T/ρ)
    f=20
    ω=2*pi*f                     #频率以及角频率
    λ=v/f                          #波长
    h=λ/25
    L=1.0                                   #绳子长度为1米
    α=0.8;c=α^2
    τ=α*h/v
    N=Int(floor(L/h))
    x=range(0,L,length=N)
    y0=zeros(N)                             #绳子的位移
    y1=zeros(N)                             #即初始时刻dy/dt=0;
    y =zeros(N)
    for k=1:200
        t=k*τ;
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=sin(ω*t)
        y0 .= y1
        y1 .= y
        plot(x,y)
        ylim([-2,2]) 
        sleep(0.1)
    end
end
rope_wave()