
# using Plots
using FileIO
using ProgressMeter

# gr(fmt=:png)
# plotly()


name_f = "50Hz_100Hz_eps9"
N_images = 14999
folder = "E:/DATA_FTP/111024/h_map_$(name_f)"
folder_save = "E:/DATA_FTP/111024/"

name = "h_map_510.jld2"
full_path = joinpath(folder,name)
h_total = load(full_path,"h_map")

Nx,Ny = size(h_total)
X = 0:1:Nx-1
Y = 0:1:Ny-1
p = 1e-2/(33)
X_ = repeat(X,1,length(Y)).*p .* 1000 # en mm
Y_ = repeat(Y',length(X),1).*p .* 1000 # en mm

""" PLOT """
# heatmap(X,Y,-h_total[:,:].*1000)
# surface(X_,Y_,-h_total[:,:].*1000, size=(800,800)) # en mm
# plot(h_total[750,10:Ny-10])

x1 = Nx÷5
x2 = Nx*2÷5
x3 = Nx*3÷5
x4 = Nx*4÷5
y1 = Ny÷5
y2 = Ny*2÷5
y3 = Ny*3÷5
y4 = Ny*4÷5

nx,ny = size(h_total)
h_profile0 = zeros(N_images,ny)
h_profile1 = zeros(N_images,ny)
h_profile2 = zeros(N_images,ny)
h_profile3 = zeros(N_images,ny)

h_profile0y = zeros(N_images,nx)
h_profile1y = zeros(N_images,nx)
h_profile2y = zeros(N_images,nx)
h_profile3y = zeros(N_images,nx)

@showprogress 1 "Progress" for i=1:N_images
    name_i = "h_map_$(i).jld2"
    full_path_i = joinpath(folder,name_i)
    h_total_i = load(full_path_i,"h_map")
    h_profile0[i,:] = h_total_i[x1,:]
    h_profile1[i,:] = h_total_i[x2,:]
    h_profile2[i,:] = h_total_i[x3,:]
    h_profile3[i,:] = h_total_i[x4,:]

    h_profile0y[i,:] = h_total_i[:,y1]
    h_profile1y[i,:] = h_total_i[:,y2]
    h_profile2y[i,:] = h_total_i[:,y3]
    h_profile3y[i,:] = h_total_i[:,y4]

end

save( joinpath(folder_save, "h_x250_$(name_f).jld2"), "h_profile", h_profile0 )
save( joinpath(folder_save, "h_x750_$(name_f).jld2"), "h_profile", h_profile1 )
save( joinpath(folder_save, "h_x1000_$(name_f).jld2"), "h_profile", h_profile2 )
save( joinpath(folder_save, "h_x1250_$(name_f).jld2"), "h_profile", h_profile3 )

save( joinpath(folder_save, "h_y300_$(name_f).jld2"), "h_profile", h_profile0y )
save( joinpath(folder_save, "h_y600_$(name_f).jld2"), "h_profile", h_profile1y )
save( joinpath(folder_save, "h_y900_$(name_f).jld2"), "h_profile", h_profile2y )
save( joinpath(folder_save, "h_y1200_$(name_f).jld2"), "h_profile", h_profile3y )

print("DONE")

#plot(h_profile[2100,:])






