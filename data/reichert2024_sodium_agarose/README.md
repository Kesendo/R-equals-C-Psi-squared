# Sodium in agarose: two figures of Reichert's thesis, read off their vector paths

**Data:** Simon Reichert's measurements, with V. Schepkin, D. Kleimaier, F. G. Zöllner and L. R. Schad; please cite his thesis, and the paper for its figures, when using these values.

**Source:** S. Reichert, *Advancing sodium triple quantum (TQ) nuclear magnetic resonance (NMR) spectroscopy and imaging*, dissertation, Heidelberg University (2024), [doi:10.11588/heidok.00034464](https://doi.org/10.11588/heidok.00034464). Two of its figures are figures of the first paper, S. Reichert, V. Schepkin, D. Kleimaier, F. G. Zöllner and L. R. Schad, *NMR in Biomedicine* 37 (2024) e5106, which prints them as raster images; the thesis draws them as vector paths, so their points can be read exactly up to the figure's own coordinates.

**Producer:** [`simulations/read_reichert_thesis_figures.py`](../../simulations/read_reichert_thesis_figures.py), run on a local copy of the thesis PDF (the PDF is not in the repository); its output is [`reader_output.txt`](reader_output.txt). Each axis is calibrated by a straight line through its tick marks; the script prints the residuals of those lines and the coordinate step of the figure's content stream, and names the curves of Fig. 4.11 by the colours of its legend.

**Used by:** [One Sodium Pool or More](../../experiments/SODIUM_ONE_POOL_OR_MORE.md), through Part 9 of [`simulations/sodium_one_pool_or_more.py`](../../simulations/sodium_one_pool_or_more.py).

| file | figure | what it holds |
|---|---|---|
| `thesis_fig3_6a_irtqtppi_fid_2pct_pc0.csv` | Fig. 3.6a, the paper's Fig. 2A | the inversion-recovery FID of the IRTQTPPI sequence along the evolution time, 2 % agarose with 154 mM NaCl, phase cycle PC0 without double-quantum suppression, 9.4 T: 640 points, time in ms, signal in the figure's normalized units |
| `thesis_fig4_11_fit_vs_points_added_2pct_pc2.csv` | Fig. 4.11, the paper's Fig. 5 | the authors' fit of their Eq. 13 (T₁s, T₁m, T₁f in ms and the T₁-TQ/SQ ratio in %) against the number of coarse points added to the dense part, 2 % agarose, PC2, 9.4 T: 13 points |

**Precision.** The time axis of Fig. 3.6a is calibrated by seven tick marks to 0.003 ms, its signal axis by eleven to 1.4·10⁻⁵; the figure's content stream writes coordinates in steps of 0.001 pt, 0.0038 ms on the time axis. The vertices lie within 0.006 ms of two evenly stepped lines, one per band, more than that step alone gives; where the rest comes from, the drawing or the time vector the authors plotted, is not known (Part 9 of the script prints it). Fig. 4.11's axes are calibrated to 0.002 points added, 1.3·10⁻⁴ ms and 2.2·10⁻⁵ %; its values carry the authors' fits, not data.
