# Import the packages
using KomaMRI, MAT, CUDA, ArgParse


"""
Script for converting processed .mat volumes into 3D Koma phantoms.
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
#print(subjects)

for subject in subjects
    # get data from input (.mat file)
    path = args["i"]*subject*"_MPM.mat"
    #print(path)
    data = MAT.matread(path)

    # define spin position arrays
    delta_x, delta_y, delta_z = data["RESOS"] * 1e-3    # resolution in mm
    x_dim, y_dim, z_dim = data["Shape"]                 # number of spins in x, y and z

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

    # get data and flatten the arrays
    x = [(x_grid...)...]
    y = [(y_grid...)...]
    z = [(z_grid...)...]
    # convert back from ms to s
    M0 = [(data["M0"]...)...]/1000
    T1 = [(data["T1"]...)...]/1000
    T2 = [(data["T2"]...)...]/1000
    T2s = [(data["T2s"]...)...]/1000
    DB = [(data["DB"]...)...]/1000

    # define the phantom
    obj = Phantom{Float64}(
            name = subject*"_MPM",
                x = x[M0.!=0],
                y = y[M0.!=0],
                z = z[M0.!=0],
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
    write_phantom(obj, args["o"]*subject*"_MPM.phantom") 
end