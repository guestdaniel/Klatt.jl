export click, plot_spectrum, calc_spectrum

# Quick click with one sample of click and 99999 samples of silence
click(n=100_000) = vcat(1.0, zeros(n-1))

# Method for plot spectrum with multiple signals
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

# Quick and dirty function to plot spectral analysis with Makie
function plot_spectrum(
    x; 
    fs=100e3, 
    xlims=(50.0, 5e3), 
    normalize=true, 
    xscale=identity, 
    annotate=false, 
    fig=Figure(),
    ax=Axis(fig[1, 1]; xscale=xscale),
)
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

    # Adjust limits and labels
    ax.xlabel = "Frequency (Hz)"
    ax.ylabel = "Power (dB)"
    xlims!(ax, xlims...)

    # Adjust things based on xscale
    if annotate
        for slope in 6.0:6.0:24.0
            x = LogRange(xlims[1] * 2, xscale == identity ? xlims[2]*0.9 : xlims[2] / 2, 500)
            y = -slope .* log2.(x ./ (xlims[1]*2))
            lines!(ax, x, y, color=:gray, linestyle=:dash)
        end
    end
    fig
end

# Quick and dirty unscaled dB power spectrum
function calc_power_spectrum(x; fs=100e3)
    f = collect(LinRange(0.0, fs, length(x)))
    X = 20 .* log10.(abs.(fft(x)))
    return f, X
end

