module Klatt

using AuditorySignalUtils
using CairoMakie
using DSP
using FFTW

export klatt_resonator, klatt_cascade, check_vw, check_resonators
export plot_spectrum, click

click(n=100_000) = vcat(1.0, zeros(n-1))

function check_vw(f0=100.0, args...; radimp=false, kwargs...)
    plot_spectrum(klatt_cascade(f0, Float64[], Float64[]; radimp=radimp, dur=5.0), args...; kwargs...)
end

function check_resonators(mode="freq")
    if mode == "freq"
        plot_spectrum(map(x -> Klatt.klatt_resonator(click(), x, 50.0), 1e3:1e3:5e3))
    else
        plot_spectrum(map(x -> Klatt.klatt_resonator(click(), 2500.0, x), [50.0, 100.0, 200.0, 400.0, 800.0]))
    end
end

function klatt_cascade(
    f0::T=100.0, 
    f::Vector{T}=[500.0, 1.1e3, 2.3e3, 4e3, 5e3], 
    bw::Vector{T}=[100.0, 200.0, 250.0, 350.0, 400.0]; 
    fs=100e3, 
    dur=1.0,
    radimp=true,
) where {T <: AbstractFloat}
    x = zeros(samples(dur, fs))
    n = samples(1/f0, fs)
    x[1:n:end] .= 1.0
    x = klatt_resonator(x, 0.0, 100.0; fs=fs)
    x = klatt_antiresonator(x, 1500.0, 6000.0; fs=fs)
    for (_f, _bw) in zip(f, bw)
        x = klatt_resonator(x, _f, _bw; fs=fs)
    end
    if radimp x = klatt_radiation_impedance(x) end
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

function klatt_antiresonator(x::Vector{T}, f::T, bw::T; fs=100e3) where {T <: AbstractFloat}
    # Calculate filter coefficients from equation ?
    C = -exp(-2π * bw/fs)
    B = 2 * exp(-π * bw/fs) * cos(2 * π * f/fs)
    A = 1 - B - C

    # Equation 4
    A′ = 1.0/A
    B′ = -B/A
    C′ = -C/A

    # Filter x
    filt([A′, B′, C′], x)
end


function klatt_radiation_impedance(x::Vector{<:AbstractFloat})
    filt([1.0, -1.0], x)
end

function plot_spectrum(
    x::Vector{<:Vector};
    xscale=identity,
    kwargs...
)
    fig=Figure()
    ax=Axis(fig[1, 1]; xscale=xscale)
    map(x) do _x
        plot_spectrum(_x; fig=fig, ax=ax, kwargs...)
    end
    fig
end

function plot_spectrum(
    x; 
    fs=100e3, 
    xlims=(50.0, 5e3), 
    normalize=true, 
    xscale=identity, 
    annotate=false, 
    mode="power",
    fig=Figure(),
    ax=Axis(fig[1, 1]; xscale=xscale),
)
    # Calculate power spectrum
    if mode == "power"
        f, X = calc_power_spectrum(x; fs=fs)
    else
        f, X = calc_magnitude_spectrum(x; fs=fs)
    end

    # Plot spectrum (only use portion within xlimits for faster plotting)
    idxs = xlims[1] .<= f .<= xlims[2]
    if normalize
        lines!(ax, f[idxs], X[idxs] .- maximum(X[idxs]))
        ylims!(ax, -80.0, 10.0)
    else
        lines!(ax, f[idxs], X[idxs])
        ylims!(ax, 20.0, 80.0)
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

function calc_magnitude_spectrum(x; fs=100e3)
    f = collect(LinRange(0.0, fs, length(x)))
    X = 10 .* log10.(abs.(fft(x)))
    return f, X
end

end # module Klatt
