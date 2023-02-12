#-----------------------------------------------------------------------------#
# RK-3格式
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

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#-----------------------------------------------------------------------------#
function numerical(nx,nt,dx,dt,x,u,α)
    un = Array{Float64}(undef, nx+1) 
    ut = Array{Float64}(undef, nx+1) 
    r = Array{Float64}(undef, nx)

    k = 1 
    

    for i = 1:nx+1
        un[i] = -sin(pi*x[i])
        u[i,k] = un[i] 
    end

    un[1] = 0.0
    un[nx+1] = 0.0

    ut[1] = 0.0
    ut[nx+1] = 0.0

    for j = 2:nt+1
        rhs(nx,dx,un,r,α)

        for i = 2:nx
            ut[i] = un[i] + dt*r[i]
        end

        rhs(nx,dx,ut,r,α)

        for i = 2:nx
            ut[i] = 0.75*un[i] + 0.25*ut[i] + 0.25*dt*r[i]
        end

        rhs(nx,dx,ut,r,α)

        for i = 2:nx
            un[i] = (1.0/3.0)*un[i] + (2.0/3.0)*ut[i] + (2.0/3.0)*dt*r[i]
        end

        k = k+1
        u[:,k] = un[:]
    end
end

function rhs(nx,dx,u,r,α)
    for i = 2:nx
        r[i] = α*(u[i+1] - 2.0*u[i] + u[i-1])/(dx*dx)
    end
end

#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
x_l = -1.0
x_r = 1.0
dx = 0.025
nx = Int64((x_r-x_l)/dx)

dt = 0.0025
t = 1.0
nt = Int64(t/dt)

α = 1/(pi*pi) # 扩散率

u = Array{Float64}(undef, nx+1, nt+1)
x = Array{Float64}(undef, nx+1)
u_e = Array{Float64}(undef, nx+1)
error = Array{Float64}(undef, nx+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1) 
    u_e[i] = -exp(-t)*sin(pi*x[i]) 
end

numerical(nx,nt,dx,dt,x,u,α)

u_error = zeros(nx+1)
u_error = u[:,nt+1] - u_e
rms_error = compute_l2norm(nx, u_error)
max_error = maximum(abs.(u_error))

# create output file for L2-norm
output = open("output2.txt", "w");
write(output, "Error details: \n");
write(output, "L-2 Norm = ", string(rms_error), " \n");
write(output, "Maximum Norm = ", string(max_error), " \n");

# create text file for final field
field_final = open("field_final2.csv", "w");
write(field_final, "x"," ", "ue", " ", "un", " ", "uerror" ," \n")

for i = 1:nx+1
    write(field_final, @sprintf("%.16f",x[i])," ",@sprintf("%.16f", u_e[i])," ",
          @sprintf("%.16f", u[i,nt+1])," ",@sprintf("%.16f", u_error[i])," \n")
end

close(field_final)
close(output);


# 绘图
using DataFrames
using CSV
using PyPlot
rc("font", family="Arial", size=16.0)


final_field = CSV.read("field_final2.csv",DataFrame)

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

fig = figure("An example", figsize=(14,6));
ax1 = fig[:add_subplot](1,2,1);
ax2 = fig[:add_subplot](1,2,2);

ax1.plot(x, u_e, lw=4, ls = "-", color="b", label="Exact solution")
ax1.plot(x, u_n, lw=4, ls = "--", color="r", label="RK3 solution")
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
#ax2.legend(fontsize=14, loc=0)

# plt[:subplot](ax1);
# plt[:subplot](ax2);

fig.tight_layout()
fig.savefig("rk3.png")
