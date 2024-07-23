using Distributed
using ProgressMeter

begin time

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

            #phase0 = DSP.unwrap( angle.( Ny0 ), dims=1:2 )
            #phase = DSP.unwrap( angle.( Ny ), dims=1:2 )

            ### test 
            #Jump = phase[Int(round(nx/2)),Int(round(ny/2))] - phase0[Int(round(nx/2)),Int(round(ny/2))]
            #if abs(Jump) >= pi
            #    N = Jump÷pi
            #    phase_unwrap = phase .- sign(Jump)*N*pi
            #else
            #    phase_unwrap = phase
            #end
            ###
            #dephase_ = phase_unwrap - phase0

            dephase_ = phase .- phase0 
            #moyenne = mean(dephase_)
            #if abs(moyenne) >= pi
            #    N = moyenne÷(pi/2)
            #    result = dephase_ .- sign(moyenne)*(N*pi/2)
            #else 
            #    result = dephase_                
            #end
            
            # cos_ = cos.(dephase_)
            # sin_ = sin.(dephase_)
            # atan.(sin_,cos_)
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

        # Confirm that lines are Sinus
        #plot(reference[250,:])
        #plot!(reference[750,:])
        #plot!(reference[250,:], seriestype=:scatter)


        resfactor_test = mean(reference)/mean(gray_)
        ref_m_gray_test = reference .- resfactor_test.*gray_

        Nx,Ny = size(ref_m_gray_test)

        # Calculate wavelength of the projected pattern
        line_ref = ref_m_gray_test[Int(floor(Nx/2)),:] # #mean(ref_m_gray_test[:,:], dims=1)
        peaks, prop = findpeaks1d(line_ref; height=0.01)


        wavelength_pix = mean( diff(peaks) )
        pspp = p*wavelength_pix

        # Multiprocessing part : pour tout k
        
        # function extract_number(filename)
        #    # Match the pattern that captures the number before the file extension
        #    MAT = match(r"(\d+)\.tiff?$", tiff_files) # # # # # # # # attention parfois .tiff ou .tif !!!
        #    print(MAT)
        #    return parse(Int, MAT.captures[1])
        # end
        
        # Sort filenames based on the extracted numeric part
        

        function parallel_analysis(start,stop,folder,folder_save,folder_def,folder_map,resfactor,gray_,ref_m_gray,th, n_, L, D, pspp)

            path_def = joinpath( folder, folder_def )
            element_in_def = readdir( path_def ) # load the name of the images

            # Import to sort with Bastler Camera
            sorted_def = Filtre(element_in_def) #sort(element_in_def, by=extract_number)

            for k = start:stop
                name_def_k = sorted_def[k]
                full_def_k = joinpath(folder,folder_def,name_def_k)
                deformed_k = channelview( FileIO.load(full_def_k) )'
                def_m_gray_k = deformed_k .- resfactor .* gray_
                DEPHASE = calculate_phase_diff_map_1D(def_m_gray_k, ref_m_gray, th, n_) 
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
        parallel_analysis(start_idx,stop_idx,folder,folder_save,folder_def,folder_map,resfactor_test,gray_,ref_m_gray_test,th, n_, L, D, pspp)
    end


    print("Running")

end
