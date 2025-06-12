# Import the packages
using KomaMRI, MAT, CUDA, ArgParse, Glob, NIfTI


"""
Script for converting processed .mat slices into 2D Koma phantoms.
"""

s = ArgParseSettings()
@add_arg_table s begin
    "-i"
        help = "Input path for processed data"
        arg_type = String
        default = "./processed_data/travelling_head"
    "-s"
        help = "Flag for using slices as phantoms"
        action = :store_true
end

args = parse_args(ARGS, s)
input_dir=args["i"]

# get path to sites
sites = filter(isdir,readdir(input_dir,join=true))
# get subject and sessions for each site
for site in sites
    # get name of site
    site = basename(site)
    @show(site)
    # get subjects
    subjects = filter(isdir,readdir(joinpath(input_dir,site),join=true)) 
    
    for subject in subjects
        # get name of subject
        subject = basename(subject)
        @show(subject)
        # get in sessions
        sessions = filter(isdir,readdir(joinpath(input_dir,site,subject),join=true))
        
        for session in sessions
            # get name of session
            session = basename(session)
            @show(session)
            input_path = joinpath(input_dir, site, subject, session)
            # get files
            files = readdir(input_path,join=true)
            
            # define target folder for phantoms to go to
            target_folder = joinpath(input_path,"phantoms_koma/")
            # if the folder does not yet exist create it
            if !isdir(target_folder)
                mkpath(target_folder)
            end

            # init all the images for error tracking
            PD_image      = nothing
            PD_file       = nothing
            T1_image      = nothing
            T2star_image  = nothing
            T2_image      = nothing
            
            # get all the processed data
            for file in files
                # load PD, T1, T2* or T2 
                if endswith(file,"PD.nii")
                    PD_file  = niread(file)
                    PD_image = niread(file).raw
                elseif endswith(file,"T1.nii")
                    T1_image = niread(file).raw
                elseif endswith(file,"T2_map.nii.gz")
                    T2_image = niread(file)#.raw
                elseif endswith(file,"T2s_OLS.nii")
                    T2star_image = niread(file).raw
                end
            end

            # check whether values are read in correctly
            @assert(!isnothing(PD_image), "No PD image found!")
            @assert(!isnothing(T1_image), "No T1 image found!")
            @assert(!isnothing(T2_image), "No T2 image found!")
            @assert(!isnothing(T2star_image), "No T2star image found!")
            @assert(size(PD_image)==size(T1_image)==size(T2_image)==size(T2star_image), "Images need to have the same size - check input images!!")

            # if the flag for slices is input
            if args["s"]
                for slice in 1:size(PD_image)[3]
                    # resolution in mm
                    delta_x = 1 #PD_file.header.pixdim[2]
                    delta_y = 1 #PD_file.header.pixdim[3]
                    #vsize = voxel_size(PD_file.header) 
                    
                    # define spin position arrays
                    x_dim = size(PD_image)[1]
                    y_dim = size(PD_image)[2]

                    # Field of Views
                    FOVx = (x_dim-1)*delta_x 
                    FOVy = (y_dim-1)*delta_y

                    # spin coordinate vectors (reshape according to shape)
                    x = (-FOVx/2):delta_x:(FOVx/2)
                    y = (-FOVy/2):delta_y:(FOVy/2) 

                    # turn into grid points
                    x_grid, y_grid = x .+ y'*0, x*0 .+ y'

                    # flatten the sliced images and set DB to zero
                    M0 = [(PD_image[1:end,1:end,slice]...)...]
                    T1 = [(T1_image[1:end,1:end,slice]...)...]
                    T2 = [(T2_image[1:end,1:end,slice]...)...]
                    T2s = [(T2star_image[1:end,1:end,slice]...)...]
                    DB = [(zeros((x_dim,y_dim))...)...]

                    # define the phantom
                    obj = Phantom{Float64}(
                            name = subject*"_"*session*"_slice"*string(slice)*"_phantom_koma",
                                x = x_grid[M0.!=0],
                                y = y_grid[M0.!=0],
                                z = 0*x_grid[M0.!=0],
                                ρ = M0[M0.!=0],
                                T1 = T1[M0.!=0],
                                T2 = T2[M0.!=0],
                                T2s = T2s[M0.!=0],
                                Δw = DB[M0.!=0], 
                                Dλ1 = DB[M0.!=0], 
                                Dλ2 = DB[M0.!=0], 
                                Dθ = DB[M0.!=0],
                        )

                    # save phantom as .phantom to output path
                    write_phantom(obj, target_folder*subject*"_"*session*"_phantom_koma_slice"*string(slice)*".phantom") 
                end    
            else
                # resolution in mm
                delta_x = 1 #PD_file.header.pixdim[2]
                delta_y = 1 #PD_file.header.pixdim[3]
                delta_z = 1 #PD_file.header.pixdim[4]
                #vsize = voxel_size(PD_file.header)   

                # define spin position arrays
                x_dim = size(PD_image)[1]
                y_dim = size(PD_image)[2]
                z_dim = size(PD_image)[3]

                # Field of Views
                FOVx = (x_dim-1)*delta_x 
                FOVy = (y_dim-1)*delta_y 
                FOVz = (z_dim-1)*delta_z 

                # spin coordinate vectors (reshape according to shape)
                x = reshape((-FOVx/2):delta_x:(FOVx/2), x_dim, 1, 1)
                y = reshape((-FOVy/2):delta_y:(FOVy/2), 1, y_dim, 1)
                z = reshape((-FOVz/2):delta_z:(FOVz/2), 1, 1, z_dim)

                # turn into grid points
                x_grid = 1*x .+ 0*y .+ 0*z
                y_grid = 0*x .+ 1*y .+ 0*z
                z_grid = 0*x .+ 0*y .+ 1*z

                # flatten positional arrays
                x = [(x_grid...)...]
                y = [(y_grid...)...]
                z = [(z_grid...)...]

                # flatten the sliced images and set DB to zero
                M0 = [(PD_image...)...]
                T1 = [(T1_image...)...]
                T2 = [(T2_image...)...]
                T2s = [(T2star_image...)...]
                DB = [(zeros((x_dim,y_dim,z_dim))...)...]

                # define the phantom
                obj = Phantom{Float64}(
                        name = subject*"_"*session*"_vol_phantom_koma",
                            x = x_grid[M0.!=0],
                            y = y_grid[M0.!=0],
                            z = 0*x_grid[M0.!=0],
                            ρ = M0[M0.!=0],
                            T1 = T1[M0.!=0],
                            T2 = T2[M0.!=0],
                            T2s = T2s[M0.!=0],
                            Δw = DB[M0.!=0], 
                            Dλ1 = DB[M0.!=0], 
                            Dλ2 = DB[M0.!=0], 
                            Dθ = DB[M0.!=0],
                    )

                # save phantom as .phantom to output path
                write_phantom(obj, target_folder*subject*"_"*session*"_phantom_koma_vol.phantom")   
            end
        end
    end
end            