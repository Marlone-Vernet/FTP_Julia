using Images
using FileIO
using ProgressMeter
using Plots
using Statistics
using DSP
using FFTW
using FindPeaks1D
plotlyjs()

function calculate_phase_diff_map_1D(Mask, Mask0, th, ns)

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
    edge_x = 20
    edge_y = 20
    nx_, ny_ = size(Mask)

    dY = Mask[edge_x:nx_-edge_x, edge_y:ny_-edge_y]
    dY0 = Mask0[edge_x:nx_-edge_x, edge_y:ny_-edge_y]
    
    nx, ny = size(dY)
    phase0 = zeros(Float64, nx, ny)
    phase = zeros(Float64, nx, ny)

    fY0 = fft( dY0, 2)
    fY = fft( dY, 2)

    imax = argmax( abs.( fY0[Int(floor(nx/2)), 9:Int(floor(ny/2))] ) )
    ifmax = imax+9
    HW = round(ifmax*th)
    W = 2*HW+1
    win = DSP.tukey(Int(W), ns)

    win2D = repeat(win', (nx))

    gauss_filt1D = zeros(Float64,nx,ny)
    gauss_filt1D[:,Int(ifmax-HW):Int(ifmax-HW+W-1)] = win2D

    # Multiplication by the filter
    Nfy0 = fY0 .* gauss_filt1D # .* for element wise Multiplication !
    Nfy = fY .* gauss_filt1D

    Ny0 = ifft(Nfy0, 2)
    Ny = ifft(Nfy, 2)

    phase0 = DSP.unwrap( DSP.unwrap( angle.( Ny0 ), dims=2), dims=1 )
    phase = DSP.unwrap( DSP.unwrap( angle.( Ny ), dims=2), dims=1 )

    dephase_ = phase .- phase0 
    return dephase_ .- mean(dephase_)
end

function Filtre(file_name) # attention filtre suivant .tif ou .tiff format a modifier
    tiff_files = filter(x -> endswith(x, ".tif"), file_name)
end

n = 100

folder = "/Users/vernet/Desktop/hydroelastic_project/test_2/"
folder_ref = "ref"
path_ref = joinpath( folder, folder_ref )
element_in_ref = Filtre( readdir( path_ref ) )

name_ref1 = element_in_ref[23]
name_ref2 = element_in_ref[24]
name_ref3 = element_in_ref[25]

full_ref1 = joinpath(folder,folder_ref,name_ref1)
reference1 = channelview( FileIO.load(full_ref1) )'

full_ref2 = joinpath(folder,folder_ref,name_ref2)
reference2 = channelview( FileIO.load(full_ref2) )'
full_ref3 = joinpath(folder,folder_ref,name_ref3)
reference3 = channelview( FileIO.load(full_ref3) )'


plot(reference1[750,:])
plot!(reference2[750,:])
plot!(reference3[750,:])


LL = length(reference1[750,:])
period1 = welch_pgram(reference1[750,:],LL)
period2 = welch_pgram(reference2[750,:],LL)
period3 = welch_pgram(reference3[750,:],LL)

psd1 = period1.power
freq1 = period1.freq * 2*pi/(0.01/40)
psd2 = period2.power
freq2 = period2.freq * 2*pi/(0.01/40)
psd3 = period3.power
freq3 = period3.freq * 2*pi/(0.01/40)

plot(freq1, psd1, xscale=:log10, yscale=:log10)
plot!(freq3, psd3, xscale=:log10, yscale=:log10)
plot!(freq2, psd2, xscale=:log10, yscale=:log10)


plot!(reference1[750,:].-reference3[750,:])

plot!(reference2[750,:])
plot!(reference3[750,:])


heatmap(reference1)
heatmap(reference2)
heatmap(reference3)


nn = 1
th = 0.5
dphase2 = calculate_phase_diff_map_1D(reference2,reference1, th, nn)
dphase3 = calculate_phase_diff_map_1D(reference3,reference1, th, nn)

heatmap(dphase2)
heatmap(dphase3)

Nx,Ny = size(reference1) # Nx : nbre de ligne , Ny : nbre de colonnes
ligne_t = zeros(n,Nx)
point_t = zeros(n)

loc_peak = zeros(n)

@showprogress 1 "Computing..." for j=1:n
    name_ref = element_in_ref[j]
    full_ref = joinpath(folder,folder_ref,name_ref)
    reference = channelview( FileIO.load(full_ref) )'

    signal = reference[:,Int(round(Ny/2))]
    ligne_t[j,:] = signal
    loc_1 = findpeaks1d(signal)[1]
    loc_peak[j] = loc_1[2]
    point_t[j] = mean(reference)
end



plot(loc_peak)


heatmap(ligne_t)

plot(ligne_t[:,750])

period = welch_pgram(point_t)

psd = period.power
freq = period.freq * 120

plot(freq, psd, xscale=:log10, yscale=:log10)



plot(point_t)


#PyPlot.figure()
#PyPlot.loglog(freq, psd)
#PyPlot.show()
