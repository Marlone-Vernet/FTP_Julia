using Plots
using FileIO
using ProgressMeter

gr(fmt=:png)
# plotly()



folder = "/home/tanu/data1/DATA_post/170624/h_20Hz_amp1/"
folder_save = "/home/tanu/data1/DATA_post/170624/" 

name = "h_map_5000.jld2"
full_path = joinpath(folder,name)
h_total = load(full_path,"h_map")

Nx,Ny = size(h_total)
X = 0:1:Nx-1
Y = 0:1:Ny-1
p = 1e-2/(40)
X_ = repeat(X,1,length(Y)).*p .* 1000 # en mm
Y_ = repeat(Y',length(X),1).*p .* 1000 # en mm


#heatmap(X,Y,-h_total[:,:].*1000)
#surface(X_,Y_,-h_total[:,:].*1000, size=(800,800)) # en mm



# MOOVIE 

function animation_film(frame)
    name_i = "h_map_$(frame).jld2"
    full_path_i = joinpath(folder,name_i)
    h_total_i = load(full_path_i,"h_map")
    surface(X_,Y_,-h_total_i.*1000,size=(800,800))
end

anim = @animate for frame in 1:600
    animation_film(frame)
end 


name_moovie = "h_t_20Hz_amp1.mp4"
path_movie = joinpath(folder_save,name_moovie)

#save(path_movie, )
