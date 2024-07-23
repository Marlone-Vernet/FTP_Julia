using Distributed
using ProgressMeter
using SharedArrays

addprocs(6)

# plotly()

# folder = "D:/Marlone/Profilometrie/150424/h_total/"
# h_total = load( joinpath(folder, "h_map_6000.jld2"), "h_map")

# Nx,Ny = size(h_total)
# Ny_ = 1250
# X = 0:1:Nx-1
# Y = 0:1:Ny_-1
# p = 1e-2/(40)
# X_ = repeat(X,1,length(Y)).*p .* 100 # en cm
# Y_ = repeat(Y',length(X),1).*p .* 100 # en cm

# surface(X_,Y_,-h_total[:,1:Ny_].*100, size=(800,800)) # en cm
# Assuming your 2D array is stored in `data`

@everywhere begin
    using SharedArrays

    using Distributed
    using Statistics
    using FileIO
    using Interpolations
    using DSP 
    using LinearAlgebra
    using FFTW
    # const SharedArray = Base.SharedArray

    function my_meshgrid(x,y) 
        xgrid = repeat(x, 1, length(y))
        ygrid = repeat(y', length(x), 1)
        return xgrid, ygrid
    end


    function spectre_2D(h_map)
        normalisation = 1 # TO BE CHANGED 
        fft2 = fft(h_map)
        fft2 = abs2.(fftshift(fft2)) / normalisation

        nx, ny = size(fft2)
        n_square = minimum((nx,ny))
        mid = div(n_square, 2)
        kx = range(-mid, mid, length=n_square)
        ky = range(-mid, mid, length=n_square)
        Kx,Ky = my_meshgrid(kx,ky)
        t = range(0, 2 * π, length=129)
        k_ = range(0, stop=kx[end], length=mid)
        k, theta = my_meshgrid(kx[mid+1:end], t[1:128])

        kxp = k .* cos.(theta)
        kyp = k .* sin.(theta)
        ii = Int(floor((nx-n_square)/2))+1:n_square + Int(floor((nx-n_square)/2))
        jj = Int(floor((ny-n_square)/2))+1:n_square + Int(floor((ny-n_square)/2))

        n1,n2 = size(kxp)
        spp = Array{Float64}(undef, n1, n2)
        
        itp = Interpolations.interpolate((kx, ky), fft2[ii,jj], Gridded(Linear()))
        # samples = hcat(vec(kxp), vec(kyp))
        for i=1:n1 
            for j=1:n2
                spp[i,j] = itp(kxp[i,j],kyp[i,j])
            end
        end

        spp_averaged = sum(spp, dims=2) .* (t[2] - t[1]) .* k_

        return collect(k_), vec(spp_averaged)
    end

    
    function process_data(folder, j)
        h_map = load(joinpath(folder, "h_map_$(j).jld2"), "h_map")
        K_, PSD_ = spectre_2D(h_map)
        # update!(p, j)  # Update progress bar
        return K_,PSD_
    end
end

# K, PSD_test = spectre_2D(h_total)
# plot(K,PSD_test, xscale=:log10, yscale=:log10, xlims=(minimum(K),maximum(K)))


function main()

    len_ = 644
    N_images = 7200
    folder_main = "D:/Marlone/Profilometrie/160424/"
    folder_h = "h_total/"
    folder = joinpath(folder_main,folder_h)
    # Parallel map over the images to compute PSD
    # p = Progress(N_images, "Processing Images", 1)
    # spectre_vs_t = SharedArray{Float64}(N_images, len_)
    #p = Progress(N_images, "Processing Images", 1)

    spectre_vs_t = pmap(j -> process_data(folder, j), 1:N_images)
    #spectre_vs_t = pmap(j -> begin
    #                            update!(p, j)  # Update progress bar
    #                            process_data(folder, j)
    #                         end, 1:N_images)
    
    #Progress.close(p)
    
    k, psd_t = unzip(spectre_vs_t)

    
    name_save = "spectre_vs_t.jld2"
    full_save = joinpath(folder_main,name_save)
    save(full_save, "k",k)
    save(full_save, "psd_t",psd_t)

    # Progress.close(p)

    return nothing
end

main()  # Execute the main function