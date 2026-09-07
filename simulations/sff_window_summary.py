"""Import-safe summaries for the raw, M^-2-normalized frequency SFF."""
import numpy as np


def sff_frequency(eigenvalues, t_arr):
    """Raw frequency SFF: |sum exp(i Im(lambda) t)|^2 / M^2."""
    mode_count = len(eigenvalues)
    freqs = eigenvalues.imag
    values = np.zeros(len(t_arr))
    for index, time in enumerate(t_arr):
        values[index] = np.abs(np.sum(np.exp(1j * freqs * time)))**2 / mode_count**2
    return values


def raw_multiset_density_scale(eigenvalues):
    """2*pi / mean adjacent gap of sorted |Im(lambda)| > 1e-10.

    Retains multiplicities: this is a multiplicity-dependent raw multiset
    density scale, not a physical spectral-resolution invariant.
    """
    frequencies = np.sort(np.abs(eigenvalues.imag[np.abs(eigenvalues.imag) > 1e-10]))
    if len(frequencies) < 2 or frequencies[-1] == frequencies[0]:
        raise ValueError("raw multiset density scale requires a nonzero frequency range")
    return 2 * np.pi / np.mean(np.diff(frequencies))


def independent_phase_reference(mode_count):
    """Diagonal independent-phase contribution, not a degeneracy-aware plateau theorem."""
    if mode_count <= 0:
        raise ValueError("mode_count must be positive")
    return 1.0 / mode_count


def summarize_windows(times, values, raw_scale):
    """Descriptive grid bins relative to the multiplicity-dependent raw scale."""
    times, values = np.asarray(times), np.asarray(values)
    below = times < 0.1 * raw_scale
    intermediate = (times > 0.1 * raw_scale) & (times < raw_scale)
    beyond = times > raw_scale
    result = {name: float(np.mean(values[mask])) if np.any(mask) else None
              for name, mask in (("below", below), ("intermediate", intermediate), ("beyond", beyond))}
    result["minimum"] = float(np.min(values[below])) if np.any(below) else None
    result["slope"] = (float(np.polyfit(times[intermediate] / raw_scale, values[intermediate], 1)[0])
                       if np.count_nonzero(intermediate) > 5 else None)
    return result


def format_sample(value):
    return "not sampled" if value is None else f"{value:.4f}"
