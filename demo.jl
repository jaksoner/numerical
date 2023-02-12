using DataFrames
using CSV
using PyPlot
rc("font", family="Arial", size=16.0)


final_field = CSV.read("field_final.csv",DataFrame)
final_field2 = CSV.read("field_final2.csv",DataFrame)
final_field3 = CSV.read("field_final3.csv",DataFrame)

x = convert(Array,final_field[:,1])
u_error = convert(Array,final_field[:,4])
u_error2 = convert(Array,final_field2[:,4])
u_error3 = convert(Array,final_field3[:,4])

for i = 1:Int64(length(u_error))
    u_error[i] = abs(u_error[i])
end
for i = 1:Int64(length(u_error2))
    u_error2[i] = abs(u_error2[i])
end
for i = 1:Int64(length(u_error3))
    u_error3[i] = abs(u_error3[i])
end

fig = figure("An example", figsize=(14,6));
ax2 = fig[:add_subplot](1,1,1);


ax2.plot(x, u_error, marker = "o", markeredgecolor="k",
        markersize=8, color="g", lw=4,label="FTCS")
ax2.plot(x, u_error2, marker = "o", markeredgecolor="k",
        markersize=8, color="r", lw=8,label="RK-3")
# ax2.plot!(x, u_error3, marker = "o", markeredgecolor="k",
#         markersize=8, color="b", lw=4,label="CN")
ax2.set_xlabel("\$x\$")
ax2.set_ylabel("\$ϵ\$")
ax2.set_title("Discretization error")
ax2.set_xlim(-1,1)
ax2.legend(fontsize=14, loc=0)

# plt[:subplot](ax1);
# plt[:subplot](ax2);

fig.savefig("error_all.png")