# Import the packages
using KomaMRI, MAT, CUDA, ArgParse, Glob


"""
Script for converting processed .mat slices into 2D Koma phantoms.
"""

s = ArgParseSettings()
@add_arg_table s begin
    "-i"
        help = "input path"
        required = true
        arg_type = String
    "-o"
        help = "output path"
        required = true
        arg_type = String
    "-s"
        help = "subject number"
        required = true
        arg_type = String
end

args = parse_args(ARGS, s)
subjects = split(args["s"],",")

for subject in subjects
    # get all the slices
    path = "sub-tle"*subject*"_MPM_Slice*.mat"
    slices = glob(path, args["i"])
    #print(size(slices))
    
    for (slice_num,slice) in enumerate(slices)
        # get data from input (.mat file)
        data = MAT.matread(slice)

        # define spin position arrays
        delta_x, delta_y, delta_z = data["RESOS"] * 1e-3    # resolution in mm
        x_dim, y_dim = data["Shape"]                 # number of spins in x, y and z

        # Field of Views
        FOVx = (x_dim-1)*delta_x 
        FOVy = (y_dim-1)*delta_y

        # spin coordinate vectors (reshape according to shape)
        x = (-FOVx/2):delta_x:(FOVx/2)
        y = (-FOVy/2):delta_y:(FOVy/2) 

        # turn into grid points
        x_grid, y_grid = x .+ y'*0, x*0 .+ y'

        # get data and flatten the arrays
        M0 = [(data["M0"]...)...]
        T1 = [(data["T1"]...)...]
        T2 = [(data["T2"]...)...]
        T2s = [(data["T2s"]...)...]
        DB = [(data["DB"]...)...]

        # define the phantom
        obj = Phantom{Float64}(
                name = subject*"_MPM",
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
        write_phantom(obj, args["o"]*subject*"_MPM_Slice"*string(slice_num-1)*".phantom") 
    end    
end