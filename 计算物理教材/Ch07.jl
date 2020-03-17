using GR

# 热传导方程（教材332）
x=0:20; a2=10;  r=a2*0.01;
u=zeros(21,25);
u[10:11,1].=1;
for j=1:24
    u[2:20,j+1]=(1-2*r)*u[2:20,j]+r*(u[1:19,j]+u[3:21,j]);
    plot(x,u[:,j]);
    ylim([0,1]);
    sleep(0.1);
end
surface(u)

# 绳波方程
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
        y[2:N-1]=copy(c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1])
        y[1]=sin(ω*t)
        y0=copy(y1)
        y1=copy(y)
        plot(x,y)
        ylim([-2,2]) 
        sleep(0.1)
    end
end