using Distributed
using ProgressMeter

@time begin

    addprocs(6)

    @everywhere using Images
    @everywhere using FileIO
    @everywhere using FFTW 
    @everywhere using DSP
    @everywhere using Statistics
    @everywhere using FindPeaks1D
    @everywhere using JLD2
    @everywhere using FixedPointNumbers
    @everywhere using FixedPointNumbers
    @everywhere using ColorTypes
    @everywhere using LinearAlgebra
    @everywhere using ProgressMeter

    @everywhere function calculate_phase_diff_map_1D(dY, dY0, th, ns)

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

        phase0 = unwrap( unwrap( angle.( Ny0 ), dims=2 ), dims=1)
        phase = unwrap( unwrap( angle.( Ny ), dims=2 ), dims=1)

        dephase_ = phase - phase0
        cos_ = cos.(dephase_)
        sin_ = sin.(dephase_)

        return atan.(sin_,cos_)
    end


    @everywhere function height_map_from_phase_map(dphase, L, D, p)
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


    @everywhere folder = "D:/Marlone/Profilometrie/150424"
    @everywhere folder_def = "deformed"
    @everywhere folder_ref = "ref"
    @everywhere folder_gray = "gray"

    @everywhere name_gray = "Basler_acA2040-120um__23597839__20240415_111412082_1.tiff"
    @everywhere name_ref = "Basler_acA2040-120um__23597839__20240415_111447911_1.tiff"

    @everywhere full_gray = joinpath(folder,folder_gray,name_gray)
    @everywhere full_ref = joinpath(folder,folder_ref,name_ref)

    @everywhere gray_ = channelview( FileIO.load(full_gray) )'
    @everywhere reference = channelview( FileIO.load(full_ref) )'

    # Confirm that lines are Sinus
    #plot(reference[250,:])
    #plot!(reference[750,:])
    #plot!(reference[250,:], seriestype=:scatter)

    @everywhere L = 1.40
    @everywhere D = 0.68
    @everywhere p = 1e-2/(40) # conversion en m/pixel
    @everywhere n_ = 1
    @everywhere th = 0.5 # taille du filtre 

    @everywhere resfactor_test = mean(reference)/mean(gray_)
    @everywhere ref_m_gray_test = reference .- resfactor_test.*gray_

    @everywhere Nx,Ny = size(ref_m_gray_test)

    # Calculate wavelength of the projected pattern
    @everywhere line_ref = ref_m_gray_test[Int(floor(Nx/2)),:] # #mean(ref_m_gray_test[:,:], dims=1)
    @everywhere peaks, prop = findpeaks1d(line_ref; height=0.01)


    @everywhere wavelength_pix = mean( diff(peaks) )
    @everywhere pspp = p*wavelength_pix

    # Multiprocessing part : pour tout k
    @everywhere folder_map = "h_total"


    @everywhere function parallel_analysis(k,folder,folder_def,folder_map,resfactor,gray_,ref_m_gray,th, n_, L, D, pspp)
        name_def_k = "Basler_acA2040-120um__23597839__20240415_111910649_$(k).tiff"
        full_def_k = joinpath(folder,folder_def,name_def_k)
        deformed_k = channelview( FileIO.load(full_def_k) )'
        def_m_gray_k = deformed_k .- resfactor .* gray_
        DEPHASE = calculate_phase_diff_map_1D(def_m_gray_k, ref_m_gray, th, n_)
        h_total = height_map_from_phase_map(DEPHASE, L, D, pspp)

        name_save = "h_map_$(k).jld2"
        full_save = joinpath(folder, folder_map, name_save)
        save(full_save, "h_map", h_total)
        return nothing
    end

    # Define the range for the loop
    @everywhere n = 21600

    print(n)

    # Parallelize the loop across workers
    @showprogress 1 "Computing..." pmap(args -> parallel_analysis(args...), [(k, folder, folder_def, folder_map, resfactor_test, gray_, ref_m_gray_test, th, n_, L, D, pspp) for k in 1:n])



    print("Job done")

end