module Tartaruga

using MAT: matread
using Statistics: mean, median
using Dates

const BASEPATH = @__DIR__

const SAMPLE_FREQ = 5000
const CALIBRATION_CONST = 470
const TIMESERIES_DOWNSAMPLE = 50

function pts_to_csv(infile)
    pts = matread(joinpath(BASEPATH, infile))["pts"]
    isarray = pts isa Array
    n = isarray ? length(pts) : length(pts["x1"])
    heave = zeros(n)
    residual = zeros(n)
    frequency = zeros(n)
    for i in Base.OneTo(n)
        if isarray
            _phi = pts[i][7, :]::Vector{Float64}
            _x1 = pts[i][1, :]::Vector{Float64}
            _x1_target = pts[i][2, :]::Vector{Float64}
        else
            _phi = vec(pts["phi_bis"][i]::Array{Float64,2})
            _x1 = vec(pts["x1"][i]::Array{Float64,2})
            _x1_target = vec(pts["x1_target"][i]::Array{Float64,2})
        end
        _x1 .*= CALIBRATION_CONST
        _x1_target .*= CALIBRATION_CONST
        # Find begining/end of each period of oscillation
        idx = findall((_phi[1:end-1] .> 6) .& (_phi[2:end] .< 0.1))
        ii = (idx[1]+1):idx[end]
        phi = _phi[ii]
        x1 = _x1[ii]
        x1_mean = mean(x1)
        x1 .-= x1_mean
        x1_target = _x1_target[ii] .- x1_mean
        # Number of periods
        n = length(idx) - 1
        # Period/frequency of oscillation
        T = mean(diff(idx))/SAMPLE_FREQ
        frequency[i] = 1/T
        # Fourier decomposition - consider the primary mode only (phase-locking removes the
        # fundamental cosine component)
        fsin = sin.(phi)
        heave[i] = sum(fsin .* x1) .* 2 / (SAMPLE_FREQ * n * T)
        residual[i] = sum(fsin .* (x1 .- x1_target)) .* 2 / (SAMPLE_FREQ * n * T)
    end
    return ["heave"=>heave, "residual"=>residual, "frequency"=>frequency]
end

function data_grid()
    return Dict(["full" => pts_to_csv("grid/19092018_V175_GridFreqAmpl_Kp200Kdneg50A0Picard_freq1.mat"),
    "stable1" => pts_to_csv("grid/19092018_V175_GridStableLCOAmpl1_Kp200Kdneg50A0Picard_freq1.mat"),
    "stable2" => pts_to_csv("grid/19092018_V175_GridStableLCOAmpl2_Kp200Kdneg50A0Picard_freq1.mat"),
    "stable3" => pts_to_csv("grid/19092018_V175_GridStableLCOAmpl3_Kp200Kdneg50A0Picard_freq1.mat"),
    "unstable1" => pts_to_csv("grid/19092018_V175_GridUnStableLCOAmpl1_Kp200Kdneg50A0Picard_freq1.mat"),
    "unstable2" => pts_to_csv("grid/19092018_V175_GridUnStableLCOAmpl2_Kp200Kdneg50A0Picard_freq1.mat"),
    "unstable3" => pts_to_csv("grid/19092018_V175_GridUnStableLCOAmpl3_Kp200Kdneg50A0Picard_freq1.mat"),])
end

function data_timeseries_uncontrolled()
    output = []
    data = matread(joinpath(BASEPATH, "20092018_V20_PhasePlot_SwitchOff_Control.mat"))
    maxlength = length(data["Heave_AfterSwitchOff_Equilibrium"])
    #
    heave = data["Heave_BeforeSwitchOff"]*10 # convert to mm
    derivative = data["DerivativeHeave_BeforeSwitchOff"]*10 # convert to mm
    time = range(0, step=1/SAMPLE_FREQ, length=maxlength)
    idx = 1:TIMESERIES_DOWNSAMPLE:length(time)
    ts = ["time"=>time[idx], "heave"=>heave[idx], "derivative"=>derivative[idx]]
    push!(output, "unstable"=>ts)
    #
    heave = data["Heave_AfterSwitchOff_StableLCO"]*10 # convert to mm
    derivative = data["DerivativeHeave_AfterSwitchOff_StableLCO"]*10 # convert to mm
    time = range(0, step=1/SAMPLE_FREQ, length=maxlength) .+ 0.054  # shift to align better (autonomous so time is arbitrary)
    idx = 1:TIMESERIES_DOWNSAMPLE:length(time)
    ts = ["time"=>time[idx], "heave"=>heave[idx], "derivative"=>derivative[idx]]
    push!(output, "lco"=>ts)
    #
    heave = data["Heave_AfterSwitchOff_Equilibrium"]*10 # convert to mm
    derivative = data["DerivativeHeave_AfterSwitchOff_Equilibrium"]*10 # convert to mm
    time = range(0, step=1/SAMPLE_FREQ, length=maxlength)
    idx = 1:TIMESERIES_DOWNSAMPLE:length(time)
    ts = ["time"=>time[idx], "heave"=>heave[idx], "derivative"=>derivative[idx]]
    push!(output, "eq"=>ts)
    return Dict(output)
end

function findamplitude(data; name="")
    ddata = diff(data)
    minima = Float64[]
    maxima = Float64[]
    # Find maxima/minima
    for i in Iterators.drop(eachindex(ddata), 1)
        if (ddata[i - 1] < 0) && (ddata[i] >= 0)
            push!(minima, data[i])
        elseif (ddata[i - 1] > 0) && (ddata[i] <= 0)
            push!(maxima, data[i])
        end
    end
    median_min = median(minima)
    median_max = median(maxima)
    amplitude = median_max - median_min
    # Discard any supurious results
    spurious_max = (maxima .> median_max + 0.2*amplitude) .|
        (maxima .< median_max - 0.2*amplitude)
    spurious_min = (minima .> median_min + 0.2*amplitude) .|
        (minima .< median_min - 0.2*amplitude)
    if any(spurious_min)
        if length(spurious_min) - count(spurious_min) < 50
            @warn "Limited data available (min)" name
        end
    end
    mean_min = mean(minima[.!spurious_min])
    if any(spurious_max)
        if length(spurious_max) - count(spurious_max) < 50
            @warn "Limited data available (max)" name
        end
    end
    mean_max = mean(maxima[.!spurious_max])
    return 0.5*(mean_max - mean_min)
end

function openloop()
    files = sort(filter(name -> contains(name, "TowardsHB") || contains(name, "BackFromHB"), readdir(joinpath(BASEPATH, "openloop"))))
    dates = Date[]
    velocities = Float64[]
    amplitudes = Float64[]
    for file in files
        meta = match(r"(........)_CBC.*_V(\d*)_", file)
        push!(dates, Date(meta[1], Dates.DateFormat("ddmmyyyy")))
        push!(velocities, parse(Int, meta[2])/10)
        push!(amplitudes, findamplitude(matread(joinpath(BASEPATH, "openloop", file))["data"][2, :]*CALIBRATION_CONST; name=file))
    end
    return (; dates, velocities, amplitudes)
end

function openloop_eq()
    files = sort(filter(name -> contains(name, "Hammer"), readdir(joinpath(BASEPATH, "openloop"))))
    velocities = Set{Float64}()
    for file in files
        # Check the first and last seconds for max/min < 0.01 to decide if an equilibrium
        heave = matread(joinpath(BASEPATH, "openloop", file))["data"][2, :]*CALIBRATION_CONST
        if (maximum(heave[1:251]) - minimum(heave[1:251])) > 5 &&
            (maximum(heave[end-250:end]) - minimum(heave[end-250:end])) > 5
            continue
        end
        meta = match(r"(........)_CBC.*_V(\d*)_", file)
        push!(velocities, parse(Int, meta[2])/10)
    end
    return velocities
end

function data_openloop()
    output = []
    lco = openloop()
    for (i, date) in enumerate(sort(unique(lco.dates)))
        idx = findall(lco.dates .== date)
        push!(output, "day$i" => ["velocity"=>lco.velocities[idx], "heave"=>lco.amplitudes[idx]])
    end
    data = sort(collect(openloop_eq()))
    push!(output, "eq" => ["velocity"=>data, "heave"=>zeros(size(data))])
    return Dict(output)
end

end
