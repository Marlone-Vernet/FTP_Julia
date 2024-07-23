using Images
using FileIO
using FFTW 
using DSP
using Statistics
using FindPeaks1D
using JLD2
using FixedPointNumbers
using ColorTypes
using LinearAlgebra
using ProgressMeter

import DSP

using Plots
plotlyjs()
#plotly()
#using Makie

#gr()#fmt=:png)



function calculate_phase_diff_map_1D(dY, dY0, th, ns)

    """
    # % Basic FTP treatment.
    # % This function takes a deformed and a reference image and calculates the phase difference map between the two.
    # %
    # % INPUTS:
    # % dY	= deformed image
    # % dY0	= reference image
    # % ns	= size of gaussian filter
    # %
    # % OUTPUT:
    # % dphase 	= phase difference map between images
    """
    nx, ny = size(dY)
    phase0 = zeros(Float64, nx, ny)
    phase = zeros(Float64, nx, ny)

    fY0 = fft( dY0, 2)
    fY = fft( dY, 2)

    imax = argmax( abs.( fY0[Int(floor(nx/2)), 9:Int(floor(ny/2))] ) )
    ifmax = imax+9
    HW = round(ifmax*th)
    W = 2*HW+1
    win = tukey(Int(W), ns)

    win2D = repeat(win', (nx))

    gauss_filt1D = zeros(Float64,nx,ny)
    gauss_filt1D[:,Int(ifmax-HW):Int(ifmax-HW+W-1)] = win2D

    # Multiplication by the filter
    Nfy0 = fY0 .* gauss_filt1D # .* for element wise Multiplication !
    Nfy = fY .* gauss_filt1D

    Ny0 = ifft(Nfy0, 2)
    Ny = ifft(Nfy, 2)

    phase0 = DSP.Unwrap.unwrap( angle.( Ny0 ), dims=1:2 )
    phase = DSP.Unwrap.unwrap( angle.( Ny ), dims=1:2 )

    dephase_ = phase - phase0

    #moyenne = mean(dephase_)
    #if moyenne >= pi
    #    N = moyenne÷(pi/2)
    #    result = dephase_ .-(N*pi/2)
    #elseif moyenne <= -pi
    #    N = moyenne÷(pi/2)
    #    result = dephase_ .+(N*pi/2) 
    #else 
    #    result = dephase_                
    #end
    # cos_ = cos.(dephase_)
    # sin_ = sin.(dephase_)

    return dephase_ .-mean(dephase_)
end


function height_map_from_phase_map(dphase, L, D, p)
    """
    Converts a phase difference map to a height map using the phase to height
    relation.
    
    INPUTS:
        dphase    = phase difference map (already unwrapped)
        L         = distance between the reference surface and the plane of the entrance  pupils
        D         = distance between centers of entrance pupils
        p         = wavelength of the projected pattern (onto the reference surface)
        spp       = physical size of the projected pixel (as seen onto the reference  surface)
        
        OUTPUT:
            h         = height map of the surface under study
"""
    return -L.*dphase./(2*pi/p*D .- dphase)
end 

#""""""""""""""""""""""""" TO BE CHANGed """"""""""""""""""""""""""""""""""
# Define the range for the loop
n = 100

folder = "/home/tanu/Bureau/FTP_DATA/030724/"
folder_def = "deformed_sweep120s_amp2"
folder_ref = "ref"
folder_gray = "gray"
folder_map = "h_test"

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

path_gray = joinpath( folder, folder_gray)
path_ref = joinpath( folder, folder_ref )

element_in_gray = readdir( path_gray )  # load the name of the images
element_in_ref = readdir( path_ref )

name_gray = element_in_gray[1]
name_ref = element_in_ref[1]

full_gray = joinpath(folder,folder_gray,name_gray)
full_ref = joinpath(folder,folder_ref,name_ref)

gray_ = channelview( FileIO.load(full_gray) )'
reference = channelview( FileIO.load(full_ref) )'

# Confirm that lines are Sinus
plot(reference[250,:])
plot!(reference[750,:])
plot!(reference[250,:], seriestype=:scatter)

#fig = Figure()
#ax = Axis(fig[1,1])
#lines(ax, reference[250,:])

L = 1.40
D = 0.68
p = 1e-2/(40) # conversion en m/pixel
n_ = 1

#list_ = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
th = 0.5
print(th)
#th = 0.5 # taille du filtre 

resfactor_test = mean(reference)/mean(gray_)
ref_m_gray_test = reference .- resfactor_test.*gray_

Nx,Ny = size(ref_m_gray_test)

# Calculate wavelength of the projected pattern
line_ref = ref_m_gray_test[Int(floor(Nx/2)),:] # #mean(ref_m_gray_test[:,:], dims=1)
peaks, prop = findpeaks1d(line_ref; height=0.01)


wavelength_pix = mean( diff(peaks) )
pspp = p*wavelength_pix


function parallel_analysis(k,folder,folder_def,folder_map,resfactor,gray_,ref_m_gray,th, n_, L, D, pspp)

    path_def = joinpath( folder, folder_def )
    element_in_def = readdir( path_def ) # load the name of the images

    name_def_k = element_in_def[k]
    full_def_k = joinpath(folder,folder_def,name_def_k)
    deformed_k = channelview( FileIO.load(full_def_k) )'
    def_m_gray_k = deformed_k .- resfactor .* gray_
    DEPHASE = calculate_phase_diff_map_1D(def_m_gray_k, ref_m_gray, th, n_) 
    h_total = height_map_from_phase_map(DEPHASE, L, D, pspp)

    name_save = "h_map_$(k).jld2"
    full_save = joinpath(folder, folder_map, name_save)
    save(full_save, "h_map", h_total)

    return DEPHASE, h_total
end

index_map = 10000
dephase, h_map = parallel_analysis(index_map,folder,folder_def,folder_map,resfactor_test,gray_,ref_m_gray_test,th, n_, L, D, pspp)




# Parallelize the loop across workers
#@showprogress 1 "Computing..." for i = 1:n

#    dephase, h_map = parallel_analysis(i,folder,folder_def,folder_map,resfactor_test,gray_,ref_m_gray_test,th, n_, L, D, pspp)
#end

Nx,Ny = size(h_map)
X = 0:1:Nx-1
Y = 0:1:Ny-1
p = 1e-2/(40)
X_ = repeat(X,1,length(Y)).*p .* 1000 # en mm
Y_ = repeat(Y',length(X),1).*p .* 1000 # en mm

plot(h_map[:,750])
#display(plot(h_map[:,750]))
#display( heatmap(X,Y,h_map[:,:].*1000))
ps = surface(X_,Y_,h_map[:,:].*1000, size=(800,800)) # en mm
plot(ps)

sigma = 2
map_smooth = imfilter(h_map, Kernel.gaussian(sigma))  # Apply Gaussian smoothing (adjust sigma as needed)

plot(h_map[:,750])
plot!(map_smooth[:,750])


ps = surface(X_,Y_,map_smooth.*1000, size=(800,800)) # en mm
plot(ps)
