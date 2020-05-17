using GR
include("parameters.jl")

# 外插方法
function absorb_interplot()
    y0=zeros(N)                             
    y1=zeros(N)                           
    y =zeros(N)
    l=y1[N-2]
    r=y1[N-1]
    for t=0:τ:32T
        y[2:N-1]=c*(y1[3:N]+y1[1:N-2])+2*(1-c)*y1[2:N-1]-y0[2:N-1]
        y[1]=gsp(t)
        y[N]=r*(2-α)-l*(1-α)
        y0 .= y1
        y1 .= y
        l=y1[N-2]
        r=y1[N-1]
        plot(X,y)
        xlim([0,L])
        ylim([-2,2])
    end
end
absorb_interplot()