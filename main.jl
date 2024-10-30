#!/usr/bin/env julia

import Pkg
try
    using JSON, QSM, NIfTI, MriResearchTools
catch
    println("[INFO] Installing Julia packages...")
    ENV["JULIA_PKG_PRECOMPILE_AUTO"] = 0
    Pkg.add(Pkg.PackageSpec(name="JSON", version=v"0.21.4"))
    Pkg.add(Pkg.PackageSpec(name="QSM", version=v"0.5.4"))
    Pkg.add(Pkg.PackageSpec(name="NIfTI", version=v"0.5.6"))
    Pkg.add(Pkg.PackageSpec(name="MriResearchTools", version=v"3.2.0"))
    using JSON, QSM, NIfTI, MriResearchTools
end

function load_json(filename)
    open(filename, "r") do file
        json_data = JSON.parse(file)
        return json_data
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

function get_B0(config_data)
    for entry in config_data["_inputs"]
        if entry["id"] == "fieldmap"
            return entry["meta"]["MagneticFieldStrength"]
        end
    end
    error("MagneticFieldStrength not found in _inputs for fieldmap")
end

function get_TE(config_data)
    for entry in config_data["_inputs"]
        if entry["id"] == "fieldmap"
            return entry["meta"]["EchoTime"]
        end
    end
    error("EchoTime not found in _inputs for fieldmap")
end

function main()
    println("[INFO] Loading config.json...")
    config_data = load_json("config.json")
    
    println("[INFO] Extracting information...")
    mask_path = config_data["mask"]
    fieldmap_path = config_data["fieldmap"]
    B0 = haskey(config_data, "B0") ? config_data["B0"] : get_B0(config_data)
    magnitude_path = haskey(config_data, "magnitude") ? config_data["magnitude"] : nothing
    magnitude_json_path = haskey(config_data, "magnitude_json") ? config_data["magnitude_json"] : nothing
    fieldmap_json_path = haskey(config_data, "fieldmap_json") ? config_data["fieldmap_json"] : nothing
    if haskey(config_data, "input_units") == false & (fieldmap_json_path != nothing)
        error("Values for either 'input_units' or 'fieldmap_json' are required!")
    end
    fieldmap_units = haskey(config_data, "input_units") ? config_data["input_units"] : load_json(fieldmap_json_path)["Units"]

    println("[INFO] Loading NIfTI images...")
    mask_nii = niread(mask_path)
    fieldmap_nii = niread(fieldmap_path)
    vsz = voxel_size_safe(fieldmap_nii.header)
    println("[INFO] Voxel size is ", vsz)

    println("[INFO] Masking fieldmap...")
    mask = !=(0).(mask_nii.raw)
    fieldmap = fieldmap_nii.raw .* mask

    println("[INFO] Converting fieldmap units to Hz...")
    γ = 267.52e6  # rad/s/T

    if fieldmap_units == "rad/s"
        @views for t in axes(fieldmap, 4)
            fieldmap[:, :, :, t] .*= inv(2 * π)
        end
    elseif fieldmap_units == "Tesla"
        @views for t in axes(fieldmap, 4)
            fieldmap[:, :, :, t] .*= γ / (2 * π)
        end
    elseif fieldmap_units == "Hz"
        println("[INFO] Fieldmap is already in Hz.")
    else
        error("Unsupported fieldmap units: $fieldmap_units")
    end

    println("[INFO] Performing V-SHARP background field removal...")
    tissue_phase, vsharp_mask = vsharp(fieldmap, mask, vsz)

    println("[INFO] Saving outputs...")
    mkpath("fmap")
    savenii(tissue_phase, "fieldmap.nii.gz", "fmap", fieldmap_nii.header)
    savenii(vsharp_mask, "mask.nii.gz", "mask", fieldmap_nii.header)

    if magnitude_path != nothing
        println("[INFO] Copying input magnitude...")
        cp(magnitude_path, joinpath("fmap", basename(magnitude_path)), force=true)
    end
    if magnitude_json_path != nothing
        println("[INFO] Copying input magnitude JSON...")
        cp(magnitude_json_path, joinpath("fmap", basename(magnitude_json_path)), force=true)
    end
    if fieldmap_json_path != nothing
        println("[INFO] Copying input fieldmap JSON...")
        cp(fieldmap_json_path, joinpath("fmap", basename(fieldmap_json_path)), force=true)
    end

    println("[INFO] V-SHARP processing completed successfully.")
end

main()

