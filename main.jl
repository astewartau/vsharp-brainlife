#!/usr/bin/env julia

import Pkg
try
    using JSON, QSM, NIfTI, MriResearchTools
catch
    println("[INFO] Installing Julia packages...")
    ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0
    Pkg.add(Pkg.PackageSpec(name="JSON", version=v"0.21.4"))
    Pkg.add(Pkg.PackageSpec(name="QSM", version=v"0.5.4"))
    Pkg.add(Pkg.PackageSpec(name="NIfTI", version=v"0.5.6"))
    Pkg.add(Pkg.PackageSpec(name="MriResearchTools", version=v"3.2.0"))
    using JSON, QSM, NIfTI, MriResearchTools
end

function load_config(filename)
	open(filename, "r") do file
		config_data = JSON.parse(file)
		return config_data
	end
end

function voxel_size_safe(header::NIfTI.NIfTI1Header)
    # Define spatial multipliers, assuming 1 (mm), 2 (meter), and 3 (micron)
    SPATIAL_UNIT_MULTIPLIERS = Dict(1 => 1.0, 2 => 1000.0, 3 => 0.001)

    # Check if `xyzt_units` contains a supported unit, default to mm if not set
    unit_key = header.xyzt_units & Int8(3)  # Bitwise AND to get unit encoding
    multiplier = get(SPATIAL_UNIT_MULTIPLIERS, unit_key, 1.0)  # Default to mm if unknown

    dim_count = min(3, header.dim[1])  # Use up to 3 spatial dimensions
    Tuple(header.pixdim[i] * multiplier for i in 2:dim_count+1)
end

function main()
    println("[INFO] Loading config.json...")
	config_data = load_config("config.json")
		
	println("[INFO] Extracting information...")
	frequency_path = config_data["frequency"]
	mask_path = config_data["mask"]
	input_units = config_data["input_units"]
	output_units = config_data["output_units"]
	B0 = config_data["B0"]
	TE = config_data["TE"]

    println("[INFO] Checking files exist...")
    all_files_exist = true
    for path in [frequency_path, mask_path]
        if !isfile(path)
            println("Error: File not found: $path")
            all_files_exist = false
        end
    end

    if !all_files_exist
        println("[ERROR] Missing files! Exiting...")
        exit(1)
    end

    println("[INFO] Loading NIfTI images...")
    frequency_nii = niread(frequency_path)
    
    vsz = voxel_size_safe(frequency_nii.header)
    println("[INFO] Read voxel size ", vsz)

    mask_nii = niread(mask_path)
    mask = !=(0).(mask_nii.raw)
    frequency = frequency_nii.raw .* mask

    println("[INFO] Converting units to Hz...")
    if input_units == "rad/s"
        frequency = frequency ./ (2π)
    elseif input_units == "radians"
        γ = 267.52e6 # rad/s/T
        @views for t in axes(frequency, 4)
            frequency[:, :, :, t] .*= inv(B0 * γ * TE)
        end
    end

    println("[INFO] Performing V-SHARP background field removal...")
    tissue_phase, vsharp_mask = vsharp(frequency, mask, vsz)

    println("[INFO] Saving outputs...")
    mkpath("tissue_freq")
    savenii(tissue_phase, "t2starw.nii.gz", "tissue_freq", frequency_nii.header)
    mkpath("mask")
    savenii(vsharp_mask, "mask.nii.gz", "mask", frequency_nii.header)

    println("[INFO] V-SHARP processing completed successfully.")
end

main()

