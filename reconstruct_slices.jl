using KomaMRI, MRIReco, NIfTI

phantom_name = "sub-tle001_MPM_Slice18"

# set target folder for raw signal and reconstructed image 
target_folder = "/home/janmeyer/KomaSimulations/Simulations/"*phantom_name*"_EPI/"
if !isdir(target_folder)
    mkpath(dirname(target_folder))
end  

# load raw signal from previous simulation (and visualize it)
f = ISMRMRDFile(target_folder*phantom_name*"_raw.h5")
raw = RawAcquisitionData(f)
#plot_signal(raw)

# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false                            # Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))               # Normalize k-space to -.5 to .5 for NUFFT
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

# reconstruct image 
image = reconstruct_2d_image(raw)
# scale to [0,1]
image = image/maximum(image)
# and save the image as .png
save_Image(target_folder*phantom_name*"_image.png", image)

# read in nifti image and replace image contents with the reconstructed simulated image
image_nii = niread("./processed_data/MPM/sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii")
niwrite(target_folder*phantom_name*"_image.nii", image_nii)
