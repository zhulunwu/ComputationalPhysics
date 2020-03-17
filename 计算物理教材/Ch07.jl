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
