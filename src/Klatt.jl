module Klatt

using AuditorySignalUtils
using CairoMakie
using DSP
using FFTW

export klatt_resonator, klatt_cascade
export plot_spectrum

function klatt_cascade(f0::T=100.0, f::Vector{T}=[500.0, 1.1e3, 2.3e3], bw::Vector{T}=[100.0, 200.0, 250.0]; fs=100e3, dur=1.0) where {T <: AbstractFloat}
    x = zeros(samples(dur, fs))
    n = samples(1/f0, fs)
    x[1:n:end] .= 1.0
    x = klatt_resonator(x, 0.0, 100.0; fs=fs)
    for (_f, _bw) in zip(f, bw)
        x = klatt_resonator(x, _f, _bw; fs=fs)
    end
    x = klatt_radiation_impedance(x)
    return x
end

function klatt_resonator(x::Vector{T}, f::T, bw::T; fs=100e3) where {T <: AbstractFloat}
    # Calculate filter coefficients from Equation 2
    C = -exp(-2π * bw/fs)
    B = 2 * exp(-π * bw/fs) * cos(2π * f/fs)
    A = 1 - B - C

    # Filter x
    filt([A], [1.0, -B, -C], x)
end

function klatt_radiation_impedance(x::Vector{<:AbstractFloat})
    filt([1.0, -1.0], x)
end

function plot_spectrum(x; fs=100e3, xlims=(50.0, 5e3), normalize=false, xscale=identity, annotate=false)
    # Create figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=xscale)

    # Calculate power spectrum
    f, X = calc_power_spectrum(x; fs=fs)

    # Plot spectrum (only use portion within xlimits for faster plotting)
    idxs = xlims[1] .<= f .<= xlims[2]
    if normalize
        lines!(ax, f[idxs], X[idxs] .- maximum(X[idxs]))
        ylims!(ax, -80.0, 10.0)
    else
        lines!(ax, f[idxs], X[idxs])
        ylims!(ax, 20.0, 80.0)
    end

    tfs = [1e3, 2e3, 4e3, 8e3]
    vals = map(tfs) do tf
        i = argmin(abs.(f .- tf))
        X[i]
    end
    map(zip(tfs, vals .- vals[1])) do (tf, val)
        println("$tf Hz | $val dB")
    end

    # Adjust limits and labels
    ax.xlabel = "Frequency (Hz)"
    ax.ylabel = "Power (dB)"
    xlims!(ax, xlims...)

    # Adjust things based on xscale
    if annotate
        for slope in 6.0:6.0:24.0
            x = LogRange(xlims[1] * 2, xlims[2] / 2, 500)
            y = -slope .* log2.(x ./ (xlims[1]*2))
            lines!(ax, x, y, color=:gray, linestyle=:dash)
        end
    end
    fig
end

function calc_power_spectrum(x; fs=100e3)
    f = collect(LinRange(0.0, fs, length(x)))
    X = 20 .* log10.(abs.(fft(x)))
    return f, X
end

end # module Klatt
