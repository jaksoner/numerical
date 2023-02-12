#-----------------------------------------------------------------------------#
# FTCS格式
#-----------------------------------------------------------------------------#
using CPUTime
using Printf
using Plots
#-----------------------------------------------------------------------------#
# 计算误差项的L-2范数
#-----------------------------------------------------------------------------#
function compute_l2norm(nx,r)
    rms = 0.0
    for i = 2:nx
        rms = rms + r[i]^2
    end
    rms = sqrt(rms/((nx-1)))
    return rms
end

x_l = -1.0 
x_r = 1.0
dx = 0.025
nx = Int64((x_r-x_l)/dx)

dt = 0.0025
t = 1.0
nt = Int64(t/dt)

α = 1/(pi*pi)

x = Array{Float64}(undef, nx+1) # 空间离散点
u_e = Array{Float64}(undef, nx+1) # 解析解
un = Array{Float64}(undef, nt+1, nx+1) # 二维数组记录每一个时刻的数值解
error = Array{Float64}(undef, nx+1) # 误差项


for i = 1:nx+1
    x[i] = x_l + dx*(i-1)  # 空间x离散点
    un[1,i] = -sin(pi*x[i]) # initial condition at t=0
    u_e[i] = -exp(-t)*sin(pi*x[i]) # 解析解
end

# 边界条件
un[1,1] = 0.0 
un[1,nx+1] = 0.0

beta = α*dt/(dx*dx)

for k = 2:nt+1
    for i = 2:nx
        un[k,i] = un[k-1,i] + beta*(un[k-1,i+1] -
                                2.0*un[k-1,i] + un[k-1,i-1])
    end
    un[k,1] = 0.0 # boundary condition at x = -1
    un[k,nx+1] = 0.0 # boundary condition at x = 1
end

# 计算误差的L-2范数
u_error = un[nt+1,:] - u_e
rms_error = compute_l2norm(nx,u_error)
max_error = maximum(abs.(u_error)) # 最大误差项的绝对值


# 保存误差相关项
output = open("output.txt", "w");
write(output, "Error details: \n");
write(output, "L-2 Norm = ", string(rms_error), " \n");
write(output, "Maximum Norm = ", string(max_error), " \n");

# 储存所有的计算量
field_final = open("field_final.csv", "w");
write(field_final, "x"," ", "ue", " ", "un", " ", "uerror" ," \n")

for i = 1:nx+1
    write(field_final, @sprintf("%.16f",x[i])," ",@sprintf("%.16f", u_e[i])," ",
          @sprintf("%.16f", un[nt+1,i])," ",@sprintf("%.16f", u_error[i])," \n")
end

close(field_final)
close(output);

# 绘图
using CSV
using PyPlot
using DataFrames
rc("font", family="Arial", size=16.0)


final_field = CSV.read("field_final.csv",DataFrame)
x = convert(Array,final_field[:,1])

u_e = convert(Array,final_field[:,2])
u_n = convert(Array,final_field[:,3])
u_error = convert(Array,final_field[:,4])

u = Array{Float64}(undef, length(u_e), 2)
u[:,1] = u_e
u[:,2] = u_n

for i = 1:Int64(length(u_error))
    u_error[i] = abs(u_error[i])
end

fig = figure("FTCS", figsize=(14,6));
ax1 = fig[:add_subplot](1,2,1);
ax2 = fig[:add_subplot](1,2,2);

ax1.plot(x, u_e, lw=4, ls = "-", color="b", label="Exact solution")
ax1.plot(x, u_n, lw=4, ls = "--", color="r", label="FTCS solution")
ax1.set_xlabel("\$x\$")
ax1.set_ylabel("\$u\$")
ax1.set_title("Solution field")
ax1.set_xlim(-1,1)
ax1.legend(fontsize=14, loc=0)

ax2.plot(x, u_error, marker = "o", markeredgecolor="k",
        markersize=8, color="g", lw=4)
ax2.set_xlabel("\$x\$")
ax2.set_ylabel("\$ϵ\$")
ax2.set_title("Discretization error")
ax2.set_xlim(-1,1)
ax2.legend(fontsize=14, loc=0)

# plt[:subplot](ax1);
# plt[:subplot](ax2);

fig.tight_layout()
fig.savefig("ftcs.png")
