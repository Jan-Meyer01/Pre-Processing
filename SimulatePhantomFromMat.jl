# Import the packages
using KomaMRI, MAT, CUDA

# create custom phantom from .mat file
file_name = "sub-tle002_MPM"
path = "/home/janmeyer/Pre-Processing/volumes_mat/MPM/"*file_name*".mat"
#file_name = "sub-tle001_MPM_Slice88.mat"
#path = "/home/janmeyer/Pre-Processing/slices_mat/MPM/"*file_name
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

"""
x_grid = [i for i in x, j in 1:size(y)[1], k in 1:size(z)[1]]
y_grid = [i for i in y, j in 1:size(x)[1], k in 1:size(z)[1]]
z_grid = [i for i in z, j in 1:size(x)[1], k in 1:size(y)[1]]


M0 = [(data["M0"]...)...]
T1 = [(data["T1"]...)...]
T2 = [(data["T2"]...)...]
T2s = [(data["T2S"]...)...]
DB = [(data["DB"]...)...]
"""
# define the phantom (flatten all arrays)
obj = Phantom{Float64}(
           name = file_name,
               x = [(x_grid...)...],
               y = [(y_grid...)...],
               z = [(z_grid...)...],
               ρ = [(data["M0"]...)...],
               T1 = [(data["T1"]...)...],
               T2 = [(data["T2"]...)...],
               T2s = [(data["T2S"]...)...],
               Δw = [(data["DB"]...)...], Dλ1 = [(data["DB"]...)...], Dλ2 = [(data["DB"]...)...], Dθ = [(data["DB"]...)...],
       )

# save phantom as .phantom to later load it into the GUI
write_phantom(obj, "/home/janmeyer/KomaSimulations/Phantoms/"*file_name*".phantom") #

# start Koma GUI and add phantom (can also load it in the GUI if it was saved before)
KomaUI()
obj_ui[] = obj

# plot phantom
#plot_phantom_map(obj, :T1)

"""
# Define scanner and sequence (phantom is already defined)
sys = Scanner()
seq = PulseDesigner.EPI_example()

# Define simulation parameters and perform simulation
sim_params = KomaMRICore.default_sim_params() 
sim_params["gpu"] = true  # use GPU for faster simulation
#raw = simulate(obj, seq, sys; sim_params)

# save raw data as a mat file
sim_params["return_type"] = "mat"
raw = simulate(obj, seq, sys; sim_params)
matwrite(joinpath("Simulated images",file_name),Dict("raw" => raw))

# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
    Nx, Ny = raw.params["reconSize"][1:2]
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    image3d  = reshape(rec.data, Nx, Ny, :)
    image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
    return image2d
end

# Perform reconstruction to get the image
image2D = reconstruct_2d_image(raw)
plot_image(image2D)

# Auxiliary function for reconstruction
function reconstruct_3d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
    Nx, Ny, Nz = raw.params["reconSize"][1:3]
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny, Nz)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    image3d  = reshape(rec.data, Nx, Ny, Nz)
    return image3d
end

# Perform reconstruction to get the image
image3D = reconstruct_3d_image(raw)
plot_image(image3D)
"""