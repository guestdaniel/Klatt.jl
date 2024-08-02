module Klatt

using AuditorySignalUtils
using CairoMakie
using DSP
using FFTW

# Handle includes
include("utilities.jl")
include("diagnostics.jl")

# Handle exports
export klatt_coef, klatt_vowel, klatt_resonator, klatt_antiresonator, klatt_cascade

# Implement core Klatt functions
"""
    klatt_coef(f, bw; fs=100e3)

Return resonator coefficients A, B, C from Equation 2 (Klatt, 1980)
"""
function klatt_coef(f, bw; fs=100e3)
    # Calculate filter coefficients from Equation 2
    C = -exp(-2π * bw/fs)
    B = 2 * exp(-π * bw/fs) * cos(2π * f/fs)
    A = 1 - B - C
    return A, B, C
end

"""
    klatt_cascade(x, f::Vector, bw::Vector; fs=100e3)

Pass time-pressure waveform `x` through cascade of resonators with specified frequencies
"""
function klatt_cascade(
    x,
    f=[500.0, 1.1e3, 2.3e3, 4e3, 5e3], 
    bw=[100.0, 200.0, 250.0, 350.0, 400.0]; 
    fs=100e3, 
)
    for (_f, _bw) in zip(f, bw)
        x = klatt_resonator(x, _f, _bw; fs=fs)
    end
    return x 
end

"""
    klatt_vowel(f0, f::Vector, bw::Vector; fs=100e3, radimp=true, dur=1.0)

Synthesize Klatt (1980) vowel with static F0 and formant frequencies.

Synthesize vowel according to simplified Klatt (1980) algorithm. First, a pulse train with
specified F0 is synthesized. This waveform is then passed through a resonator with f=0 Hz
and bw=100 Hz, an antiresonator with f=1500 Hz and bw=6000 Hz, a cascade of resonators with
the specified frequencies and bandwidths, and (optionally) a filter representing acoustic
radiation impedance at the lips.
"""
function klatt_vowel(
    f0=100.0, 
    f=[500.0, 1.1e3, 2.3e3, 4e3, 5e3], 
    bw=[100.0, 200.0, 200.0, 350.0, 400.0]; 
    fs=100e3, 
    radimp=true, 
    dur=1.0,
)
    x = zeros(samples(dur, fs))
    n = samples(1/f0, fs)
    x[1:n:end] .= 1.0
    x = klatt_resonator(x, 0.0, 100.0; fs=fs)
    x = klatt_antiresonator(x, 1500.0, 6000.0; fs=fs)
    x = klatt_cascade(x, f, bw; fs=fs)
    if radimp x = klatt_radiation_impedance(x) end
    return x
end

"""
    klatt_resonator(x, f, bw; fs=100e3)

Applies resonator with frequency `f` and bandwidth `bw` to signal `x`
"""
function klatt_resonator(x, f, bw; fs=100e3)
    A, B, C = klatt_coef(f, bw; fs=fs)
    filt([A], [1.0, -B, -C], x)
end

"""
    klatt_antiresonator(x, f, bw; fs=100e3)

Applies antiresonator with frequency `f` and bandwidth `bw` to signal `x`
"""
function klatt_antiresonator(x, f, bw; fs=100e3)
    A, B, C = klatt_coef(f, bw; fs=fs)
    filt([1.0/A, -B/A, -C/A], x)
end

function klatt_radiation_impedance(x)
    filt([1.0, -1.0], x)
end

end # module Klatt
