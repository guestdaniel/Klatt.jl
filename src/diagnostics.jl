export check_vw, check_resonators, check_antiresonators

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

function check_antiresonators(mode="freq")
    if mode == "freq"
        plot_spectrum(map(x -> Klatt.klatt_antiresonator(click(), x, 50.0), 1e3:1e3:5e3))
    else
        plot_spectrum(map(x -> Klatt.klatt_antiresonator(click(), 2500.0, x), [50.0, 100.0, 200.0, 400.0, 800.0]))
    end
end

