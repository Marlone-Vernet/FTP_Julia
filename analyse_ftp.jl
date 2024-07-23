using Distributed
using ProgressMeter

@time begin

    addprocs(12)

    @everywhere begin
        using Images
        using FileIO
        using FFTW 
        import DSP
        using Statistics
        using FindPeaks1D
        using JLD2
        using FixedPointNumbers
        using ColorTypes
        using LinearAlgebra
        using ProgressMeter

        edge_x = 20
        edge_y = 20

        function calculate_phase_ref(Mask0, th, ns)
            " This function compute the phase map of the reference and is used in calculate_phase_diff_map_1D "
            nx_, ny_ = size(Mask0)
            reference_mask = Mask0[edge_x:nx_-edge_x, edge_y:ny_-edge_y]

            nx, ny = size(reference_mask)
            phase0 = zeros(Float64, nx, ny)
            fft2D_reference = fft( reference_mask, 2)

            """ compute the filter """
            imax = argmax( abs.( fft2D_reference[Int(floor(nx/2)), 9:Int(floor(ny/2))] ) )
            ifmax = imax+9
            HW = round(ifmax*th)
            W = 2*HW+1
            win = DSP.tukey(Int(W), ns)
            win2D = repeat(win', (nx))
            gauss_filt1D = zeros(Float64,nx,ny)
            gauss_filt1D[:,Int(ifmax-HW):Int(ifmax-HW+W-1)] = win2D

            """ apply filter to reference """
            fft2D_reference_filtered = fft2D_reference .* gauss_filt1D # .* for element wise Multiplication !
            real_ref = ifft(fft2D_reference_filtered, 2) # return in real space to compute phase map
            phase_ref = DSP.unwrap( DSP.unwrap( angle.( real_ref ), dims=2), dims=1 )

            return gauss_filt1D, phase_ref
        end 

        function calculate_phase_diff_map_1D(Mask, gaussian_filter, phase_ref)

            """
            # % Basic FTP treatment.
            # % This function takes a deformed and a reference image and calculates the phase difference map between the two.
            # %
            # % INPUTS:
            # % Mask = deformed image
            # % gaussian_filter = gaussian filter, computed in calculate_phase_ref
            # % phase_ref = map of the reference image phase, computed in calculate_phase_ref
            # %
            # % OUTPUT:
            # % dphase_ 	= phase difference map between images
            """

            nx_, ny_ = size(Mask)
            deformed_mask = Mask[edge_x:nx_-edge_x, edge_y:ny_-edge_y]

            nx, ny = size(deformed_mask)
            phase = zeros(Float64, nx, ny)

            fft_2D_deformed = fft( deformed_mask, 2)

            # Multiplication by the filter
            fft_2D_deformed_filter = fft_2D_deformed .* gaussian_filter
            real_deformed = ifft(fft_2D_deformed_filter, 2)

            phase = DSP.unwrap( DSP.unwrap( angle.( real_deformed ), dims=2), dims=1 )
            dephase_ = phase .- phase_ref 

            return dephase_ .- mean(dephase_)
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

        function Filtre(file_name) # attention filtre suivant .tif ou .tiff format a modifier
            tiff_files = filter(x -> endswith(x, ".tif"), file_name)
        end

        # TO BE CHECKED !!
        p = 1e-2/(37) # conversion en m/pixel
        L = 1.46
        D = 0.64
        n_ = 1
        th = 0.5 # taille du filtre 
	    sigma = 2 # Gaussian filter
        # Define the range for the loop
        n = 100

        folder = "/Users/vernet/Desktop/hydroelastic_project/test_2/" #"/home/tanu/Bureau/FTP_DATA/180724/"
        folder_save =  "/Users/vernet/Desktop/hydroelastic_project/test_2/" #"/home/tanu/data1/DATA_post/180724/"
        folder_def = "def"
        folder_ref = "ref"
        folder_gray = "gray"
        folder_map = "h_map"

        """ INITIALISATION """
        path_gray = joinpath( folder, folder_gray )
        path_ref = joinpath( folder, folder_ref )

        element_in_gray = Filtre( readdir( path_gray ) )  # load the name of the images
        element_in_ref = Filtre( readdir( path_ref ) )

        name_gray = element_in_gray[1]
        name_ref = element_in_ref[1]

        full_gray = joinpath(folder,folder_gray,name_gray)
        full_ref = joinpath(folder,folder_ref,name_ref)

        gray_ = channelview( FileIO.load(full_gray) )'
        reference = channelview( FileIO.load(full_ref) )'

        resfactor_test = mean(reference)/mean(gray_)
        ref_m_gray_test = reference .- resfactor_test.*gray_

        Nx,Ny = size(ref_m_gray_test)

        # Calculate wavelength of the projected pattern
        line_ref = ref_m_gray_test[Int(floor(Nx/2)),:] # #mean(ref_m_gray_test[:,:], dims=1)
        peaks, prop = findpeaks1d(line_ref; height=0.01)
        wavelength_pix = mean( diff(peaks) )
        pspp = p*wavelength_pix

        """ Compute phase_ref & gaussian_filter - both used in  """
        gaussian_filter, phase_ref = calculate_phase_ref(ref_m_gray_test, th, n_)

        # function extract_number(filename)
        #    # Match the pattern that captures the number before the file extension
        #    MAT = match(r"(\d+)\.tiff?$", tiff_files) # # # # # # # # attention parfois .tiff ou .tif !!!
        #    print(MAT)
        #    return parse(Int, MAT.captures[1])
        # end
        # Sort filenames based on the extracted numeric part
        
        # Multiprocessing part : pour tout k

        function parallel_analysis(start,stop,folder,folder_save,folder_def,folder_map,resfactor,gray_,gaussian_filter,phase_ref, L, D, pspp)

            path_def = joinpath( folder, folder_def )
            element_in_def = readdir( path_def ) # load the name of the images

            # Import to sort with Bastler Camera
            sorted_def = Filtre(element_in_def) #sort(element_in_def, by=extract_number)

            for k = start:stop
                name_def_k = sorted_def[k]
                full_def_k = joinpath(folder,folder_def,name_def_k)
                deformed_k = channelview( FileIO.load(full_def_k) )'
                def_m_gray_k = deformed_k .- resfactor .* gray_
                DEPHASE = calculate_phase_diff_map_1D(def_m_gray_k, gaussian_filter, phase_ref) 
                h_total = height_map_from_phase_map(DEPHASE, L, D, pspp)

                h_total_filter = imfilter(h_total, Kernel.gaussian(sigma))  # Apply Gaussian smoothing (adjust sigma as needed)
                
                name_save = "h_map_$(k).jld2"
                full_save = joinpath(folder_save, folder_map, name_save)
                save(full_save, "h_map", h_total_filter)
            end
            return nothing
        end
    end 

    chunk_size = n ÷ nworkers()
    remainder = n % nworkers()

    print(n)

	N_python = n-1
    # Parallelize the loop across workers
    @distributed for i = 1:chunk_size:N_python 
        start_idx = i
        stop_idx = min(i + chunk_size - 1, N_python)
        if remainder > 0 && i == N_python - remainder + 1
            stop_idx = N_python  # Assign remaining iterations to the last worker
        end
        parallel_analysis(start_idx,
                        stop_idx,
                        folder,
                        folder_save,
                        folder_def,
                        folder_map,
                        resfactor_test,
                        gray_,
                        gaussian_filter,
                        phase_ref,
                        L, D, pspp)
    end


    print("Running")

end


print("Done")