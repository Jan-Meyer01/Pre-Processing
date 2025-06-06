using KomaMRI, CUDA, MRIReco, NIfTI

# Script to simulate multiple slices and store them into a single 3D nifti image afterwards

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

# Define scanner and sequence
sys = Scanner()
seq = read_seq("/home/janmeyer/KomaSimulations/Sequences/EPI_SE.seq")
#seq = read_seq("/home/janmeyer/KomaSimulations/Sequences/EPI_SE.seq")

# Define simulation parameters and add the GPU device (number 0 on the server)
sim_params = KomaMRICore.default_sim_params() 
sim_params["gpu"] = true
sim_params["gpu_device"] = 0

# define phantom name and slice numbers
phantom_name = "sub-tle001_MPM"
slice_num = range(start=20,stop=155,step=1)

# set target folder for raw signal and reconstructed image 
target_folder = "/home/janmeyer/KomaSimulations/Simulations/"*phantom_name*"_EPI_ComparisonJEMRIS/"
if !isdir(target_folder)
    mkpath(dirname(target_folder))
end  

# init array for images
images = []

for i in slice_num
    obj = read_phantom("/home/janmeyer/KomaSimulations/Phantoms/Slices/"*phantom_name*"_Slice"*string(i)*".phantom")

    # start the simulation and save raw signal as .h5
    raw = simulate(obj, seq, sys; sim_params)
    file = ISMRMRDFile(target_folder*phantom_name*"_Slice"*string(i)*"_raw.h5")
    save(file, raw)

    # Perform reconstruction and save the image as .png
    push!(images,reconstruct_2d_image(raw))
    # scale to [0,1]
    #image = image/maximum(image)
    # and save the image as .png
    #save_Image(target_folder*phantom_name*"_image.png", image)
end

# read in nifti image and replace image contents with the reconstructed simulated image
#subject_num = split(phantom_name,"_")
#image_nii = niread("./processed_data/MPM/"*subject_num[1]*"_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii")
volume = zeros(size(images[1],1),size(images[1],2),length(slice_num))
for i in 1:length(slice_num)
    volume[:,:,i] = images[i]
end    
#image_nii = volume
image_nii = NIVolume(volume)
niwrite(target_folder*phantom_name*"_imageKoma.nii", image_nii)

# additionally export the scanner, sequence, raw signal and image data to .mat files
#export_2_mat(seq, obj, sys, file, rec_params, image, target_folder; type="custom", matfilename=phantom_name)



# functions for exporting to Matlab 
function export_2_mat_sequence(seq, matfolder; matfilename="seq_sequence.mat")

    # Get the samples of the events in the sequence
    samples = get_samples(seq)

    # Create a dictionary with the event samples of the sequence
    Gx = hcat(samples.gx.t, samples.gx.A)
    Gy = hcat(samples.gy.t, samples.gy.A)
    Gz = hcat(samples.gz.t, samples.gz.A)
    RF_AM = hcat(samples.rf.t, samples.rf.A)
    RF_FM = hcat(samples.Δf.t, samples.Δf.A)
    ADCS = hcat(samples.adc.t, samples.adc.A)
    seq_dict = Dict(
        "Gx" => Gx,
        "Gy" => Gy,
        "Gz" => Gz,
        "RF_AM" => RF_AM,
        "RF_FM" => RF_FM,
        "ADCS" => ADCS
    )

    # Write to matlab file
    matwrite(joinpath(matfolder, matfilename), Dict("sequence" => seq_dict))
end

function export_2_mat_kspace(seq, matfolder; matfilename="seq_kspace.mat")
    kspace, kspace_adc = get_kspace(seq; Δt=1)
    matwrite(joinpath(matfolder, matfilename), Dict("kspace" => kspace, "kspace_adc" => kspace_adc))
end

function export_2_mat_moments(seq, matfolder; matfilename="seq_moments.mat")
    t, Δt = KomaMRIBase.get_variable_times(seq; Δt=1)
    t = t[1:end-1]
    k0, _ =  KomaMRIBase.get_kspace(seq; Δt=1)
    k1, _ =  KomaMRIBase.get_M1(seq; Δt=1)
    k2, _ =  KomaMRIBase.get_M2(seq; Δt=1)
    moments = hcat(t, k0, k1, k2)
    matwrite(joinpath(matfolder, matfilename), Dict("moments" => moments))
end

function export_2_mat_phantom(phantom, matfolder; matfilename="phantom.mat")
    phantom_dict = Dict("name" => phantom.name,
                "columns" => ["x", "y", "z", "rho", "T1", "T2", "T2s", "delta_omega"],
                "data" => hcat(phantom.x, phantom.y, phantom.z, phantom.ρ, phantom.T1, phantom.T2, phantom.T2s, phantom.Δw))
    matwrite(joinpath(matfolder, matfilename), Dict("phantom" => phantom_dict))
end

function export_2_mat_scanner(sys, matfolder; matfilename="scanner.mat")
    sys_dict = Dict("B0" => sys.B0,
                "B1" => sys.B1,
                "Gmax" => sys.Gmax,
                "Smax" => sys.Smax,
                "ADC_dt" => sys.ADC_Δt,
                "seq_dt" => sys.seq_Δt,
                "GR_dt" => sys.GR_Δt,
                "RF_dt" => sys.RF_Δt,
                "RF_ring_down_T" => sys.RF_ring_down_T,
                "RF_dead_time_T" => sys.RF_dead_time_T,
                "ADC_dead_time_T" => sys.ADC_dead_time_T)
    matwrite(joinpath(matfolder, matfilename), Dict("scanner" => sys_dict))
end

function export_2_mat_raw(raw_ismrmrd, matfolder; matfilename="raw.mat");
    if haskey(raw_ismrmrd.params, "userParameters")
        dict_for_mat = Dict()
        dict_user_params = raw_ismrmrd.params["userParameters"]
        for (key, value) in dict_user_params
            if key == "Δt_rf"
                dict_for_mat["dt_rf"] = value
            elseif key == "Δt"
                dict_for_mat["dt_rf"] = value
            else
                dict_for_mat[key] = value
            end
        end
        matwrite(joinpath(matfolder, "sim_params.mat"), dict_for_mat)

        not_Koma = raw_ismrmrd.params["systemVendor"] != "KomaMRI.jl"
        t = Float64[]
        signal = ComplexF64[]
        current_t0 = 0
        for p in raw_ismrmrd.profiles
        	dt = p.head.sample_time_us != 0 ? p.head.sample_time_us * 1e-3 : 1
        	t0 = p.head.acquisition_time_stamp * 1e-3 #This parameter is used in Koma to store the time offset
            N =  p.head.number_of_samples != 0 ? p.head.number_of_samples : 1
            if not_Koma
        		t0 = current_t0 * dt
                current_t0 += N
            end
            if N != 1
                append!(t, t0.+(0:dt:dt*(N-1)))
            else
                append!(t, t0)
            end
            append!(signal, p.data[:,1]) #Just one coil
            #To generate gap
            append!(t, t[end])
            append!(signal, [Inf + Inf*1im])
        end
        raw_dict = hcat(t, signal)
        matwrite(joinpath(matfolder, matfilename), Dict("raw" => raw_dict))
    end

end

function export_2_mat_image(image, rec_params, matfolder; matfilename="image.mat")
    if haskey(rec_params, :reconSize)
        dict_rec_params = Dict("reco" => rec_params[:reco],
                            "Nx" => rec_params[:reconSize][1],
                            "Ny" => rec_params[:reconSize][2])
        matwrite(joinpath(matfolder, "rec_params.mat"), Dict("rec_params" => dict_rec_params))
    end

    matwrite(joinpath(matfolder, matfilename), Dict("image" => image))
end

function export_2_mat(seq, phantom, sys, raw_ismrmrd, rec_params, image, matfolder; type="all", matfilename="data.mat")
    head = matfilename
    if type=="all"
        export_2_mat_sequence(seq, matfolder)
        export_2_mat_kspace(seq, matfolder)
        export_2_mat_moments(seq, matfolder)
        export_2_mat_phantom(phantom, matfolder)
        export_2_mat_scanner(sys, matfolder)
        export_2_mat_raw(raw_ismrmrd, matfolder)
        export_2_mat_image(image, rec_params, matfolder)
    elseif type=="sequence"
		export_2_mat_sequence(seq, matfolder; matfilename=(head*"_sequence.mat"))
        export_2_mat_kspace(seq, matfolder; matfilename=(head*"_kspace.mat"))
        export_2_mat_moments(seq, matfolder; matfilename=(head*"_moments.mat"))
	elseif type=="phantom"
		export_2_mat_phantom(phantom, matfolder; matfilename=(head*"_phantom.mat"))
    elseif type=="scanner"
		export_2_mat_scanner(sys, matfolder; matfilename=(head*"_scanner.mat"))
    elseif type=="raw"
		export_2_mat_raw(raw_ismrmrd, matfolder; matfilename=(head*"_raw.mat"))
    elseif type=="image"
		export_2_mat_image(image, rec_params, matfolder; matfilename=(head*"_image.mat"))
    elseif type=="custom"
		export_2_mat_sequence(seq, matfolder; matfilename=(head*"_sequence.mat"))
		export_2_mat_scanner(sys, matfolder; matfilename=(head*"_scanner.mat"))
		export_2_mat_raw(raw_ismrmrd, matfolder; matfilename=(head*"_raw.mat"))
		export_2_mat_image(image, rec_params, matfolder; matfilename=(head*"_image.mat"))
	end
end