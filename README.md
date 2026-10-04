# EPG-FLAIR

Adaptation of TSE EPG code to include FLAIR, plus a FLAIR EPI test workflow. This is a continuation of [EPG-X](https://github.com/mriphysics/EPG-X?utm_source=chatgpt.com), extending TSE simulation to TSE FLAIR and then to EPI FLAIR using the EPG structure and principles learned from [EPG repository](https://github.com/imr-framework/epg?utm_source=chatgpt.com).

This version is an adaption of my KURF partner's work Valeriia Ognevaia (https://github.com/viilerr/3D-EPG-EPI-FLAIR), looking for any gaps and concepts which can be added to make the code theoretically more robust. This code could be improved further.

## Summary


## Segmentation

**The old `Test_EPG_FLAIR_Explained.m` simulated ONE 64-echo EPI train, let it decay all the way down, then copied that same decayed data 3 times to show there were 3 segments (`fill_epi2d_kspace` just pasted the same array I believe into each block). 

Real segmentation is supposed to reset the signal with a fresh excitation for each segment — I believe the old code didn't fully do that, so Figures 2 and 5 could not show the actual benefit of segmenting an EPI train.

**What I changed**

- A new nSegments option. simulate_epi3d now loops over segments; each segment reruns its own FLAIR prep and its own excitation from full magnetization, then only decays over its own short share of the echoes,
instead of one long decay curve covering everything. New outputs were added to info: nSegments, segLengths, segmentIndex, echoTimeInSegment, and
per-segment versions of MzAfterInversion/MzBeforeExcitation to correct it.
- Figure 2 now plots the real per-segment reset (dotted lines mark each new excitation). Figure 5's k-space colouring now uses each segment's true simulated signal instead of one tiled copy. I also removed the old fill_epi2d_kspace tiling function as I didn't believe it was necessary nor I could figure its use properly


## Partition

The partitions still used the pre-existing shortcut. Fn(kyOrder,kzOrder(iz)) = F0(:)*partitionDecay(iz). This computed one ky signal train once and pasted it into every partition,
scaled only by a single scalar per partition (`1` for every partition by default, since partitionSignalDecay defaults to false). From this it led to no fresh inversion,excitation,and partition-specific decay.

**What I changed**

** I made it so partitions now get their own independent FLAIR preparation and excitation, exactly like segments already do. In `EPG_FLAIR.m` / `EPG_FLAIR_Explained.m`, a new outer `for iz = 1:nPartitions`
loop wraps the existing segment loop: for a given partition, its nPE ky lines are still filled across nSegments shots (unchanged), and once that partition's ky-plane is complete, iz increments and the whole process reruns from a fresh excitation for the next partition. 


## Timing estimates

Based on info.shotDurationMs (per segment: T2prepTE + TI + that segment's own readout time) and
info.totalScanTimeMs (nPartitions * sum(shotDurationMs)), I made it so Test_EPG_FLAIR_Explained.m prints a time estimate to the command window.

With the default test-script numbers this comes out to roughly **5 minutes** for 96 shots (`3` segments `x` `32` partitions)


## An attempt at EPI vs TSE

`EPI_vs_TSE_comparison.m` was made to compare readout strategies in terms of FLAIR contrast efficiency by running the same flair prep through 2 different readouts for each tissue.

- **TSE (refocused):** each echo is refocused by a 180 degree CPMG pulse, signal decays with `T2` across the echo train (slow decay rate).

- **EPI (gradient echo):** nothing refocuses off-resonance dephasing between echoes, so signal decays with `T2star`

Both readouts share the identical flair prep code path, any difference in the results comes from the readout physics itself, not a mismatched preparation.

A single bar chart comparing WM-GM contrast and CSF nulling at the effective echo, TSE vs. EPI, both at a shared `ETL=64`.

However, TSE actually loses *more* total signal by the last echo than EPI does, despite `T2` decay being much slower per millisecond than `T2star` decay. I believe it may be due to the TSE's realistic echo spacing is longer than EPI's, so the same 64-echo train may takes TSE far longer in real time to acquire — long enough that the slower *rate* of decay is outweighed by the much longer *duration* of decay. TSE's slower decay rate does not guarantee less signal by itself.


**Tissue parameters corrected against 7T online literature:** I also changed the placeholder values for WM and GM while building this comparison, the WM and GM `T1`/`T2` values used throughout the project were found to be slightly different. When I searched online for these values, I tried to correct these 7T values to reflect from what I found in the literature:


## New Files

- `EPI_vs_TSE_comparison.m` directly compares the EPI and TSE readouts for
FLAIR contrast efficiency — the core project deliverable, see "Edit 4"
above.
- `FUTURE_WORK_STEADY_STATE.md`. - to be explained later.

## 3D FLAIR EPI model


## Figure guide

`Test_EPG_FLAIR_Explained.m` uses `T2prepTE = 50 ms`, so the recovery curves do not start at
`-1`. The T2-prep block attenuates longitudinal magnetisation by
`exp(-T2prepTE/T2)` before inversion, then inversion flips it negative. With
the current tissue values:

- WM starts at `Mz = -exp(-50/39) = -0.277`.
- GM starts at `Mz = -exp(-50/55) = -0.403`.
- CSF starts at `Mz = -exp(-50/1200) = -0.959`.

These values are all above `-1`, and Figure 1 marks the starting points plus
a `-1` reference line so the T2-prep effect is visible.

- Figure 1, `FLAIR + EPI`: one continuous timeline for segment 1.

- Figure 2, `EPI readout`: the full `nPE`-line echo train (all `nSegments`
shots concatenated) against echo number/ky line. Dotted vertical lines mark
each new segment's excitation, so the T2*-driven decay resets and restarts
at every segment boundary instead of decaying monotonically across the
whole train.

- Figure 3, `Dynamic EPG states`: tracks the longitudinal `Z0` state during
the EPI readout after FLAIR preparation. The same per-segment resets from
Figure 2 are also visible here.

- Figure 4, `Effective echo`: compares WM, GM, and CSF signal at the k-space
centre echo.

- Figure 5, `2D k-space`: fills a 2D zigzag EPI trajectory for WM and GM using
`nSegments x linesPerSegment` lines (default `3 x 64 = 192`), coloured by the
true simulated signal for the segment that acquired each line.

- Figure 6, `PSF`: computes the 2D point-spread function from the filled 2D
k-space to show the spatial blurring caused by EPI signal weighting.

- Figure 7, `PSF k-space`: shows the 2D PSF in k-space before the image-domain
point spread function is calculated.

- Figure 8, `TI optimisation`: sweeps TI to compare WM-GM contrast and CSF
suppression.


Below are the eight figures produced by the simulation. These are updated from the original repo from Valeriia Ognevaia, please refer to her code for better explanations.

## Figure 1: FLAIR Recovery with EPI

<img width="1145" height="784" alt="Image" src="https://github.com/user-attachments/assets/41f92d0d-367e-4459-8c90-1a6cf945cf6f" />

This plot charts the T2-prepared FLAIR inversion recovery curve in the EPI sequence as longitudinal magnetisation recovers after inversion. The sequence includes the EPI readout after excitation, so the recovery curve is followed by the echo train response rather than a purely TSE-only timing model. The curve should remain close to the analytic form `Mz(t) = 1 - 2*exp(-t/T1)` before excitation, and the null point should still occur near `TI ≈ T1 * ln(2)` in the simple case. This figure validates that the code applies the inversion pulse correctly, computes `Mz_init` correctly, and then feeds that recovered `Z0` state into the subsequent EPI simulation.

## Figure 2: EPI Echo Train

<img width="1145" height="784" alt="Image" src="https://github.com/user-attachments/assets/3183bad6-6dca-47ec-8544-ca96eb0d8811" />

This figure plots the simulated echo magnitudes measured during the EPI readout across consecutive echoes, now spanning all `nSegments` shots concatenated. The expected signature is a sawtooth pattern: each segment jumps back up to a fresh excitation amplitude and then decays by `T2star` (not `T2` — this is a gradient-echo EPI readout, so there is no 180° refocusing pulse mid-train) over its own shorter echo count, rather than one long monotonic decay across the entire train. To verify the simulation, check that echo spacing within a segment matches `ESP`, that amplitudes decay as predicted for the chosen `T2star`, that each segment resets at its own excitation, and that phase progression matches the RF phase schedule. If diffusion or off resonance is enabled, their effects (attenuation or phase shifts) will also appear in the echo train. Agreement with these expectations demonstrates the simulation correctly applies evolution between RF events, restarts correctly at each segment boundary, and records echoes at the intended times.

**Note on realism:** with `nSegments=1` (a full single-shot train), decay is governed purely by `T2star`, which is much shorter than `T2` (no refocusing). At clinically realistic resolutions this decays to near-zero signal well before the last echo — that is expected, correct physics for an unsegmented gradient-echo EPI train, not a simulation bug. It demonstrates why real 3D-EPI FLAIR acquisitions use segmentation (and typically also parallel imaging/partial Fourier, not yet modelled here) rather than one long single-shot train.

## Figure 3: Dynamic EPG States

<img width="1145" height="784" alt="Image" src="https://github.com/user-attachments/assets/68c2b813-5c1c-4745-b3de-5a8e02226f0e" />

This plot shows how coherence orders (transverse `F` states and longitudinal `Z` states) evolve through the pulse train. The horizontal axis represents pulse or time index and the vertical axis labels EPG order; intensity encodes magnitude. Important validation points are: (1) higher order transverse coherences appear immediately after refocusing pulses and then decay through T2; (2) the `Z0` longitudinal state displays inversion and recovery behaviour when `TI` is applied (for FLAIR); and (3) conjugation relationships between positive and negative coherence orders are preserved. Observing these behaviours confirms that the `RF_rot` rotation, the replication `build_T`, the `S` shift matrices, and the composite relax shift operator `SE` are operating together correctly to evolve the state vector.

## Figure 4: Effective echo

<img width="1145" height="784" alt="Image" src="https://github.com/user-attachments/assets/aa6e405f-5db9-4aec-bc62-03dfb03883d8" />

This figure compares the effective echo amplitude for WM, GM, and CSF at the central k-space echo of the EPI trajectory. It shows how the signal weighting at the echo centre changes with tissue specific `T1`, `T2`, and `T2*` values and how that translates into the contrast seen in the reconstructed image. Inspect the relative peak amplitudes, the central decay behaviour, and the ordering of the curves across tissues to confirm that the EPI readout is being weighted correctly before the k-space data are assembled.

## Figure 5: 2D WM/GM k-space

<img width="1315" height="867" alt="Image" src="https://github.com/user-attachments/assets/01ce383b-ff0f-415d-8412-5902d3bbd7df" />

Figure 5 shows the 2D EPI k-space trajectory for WM and GM separately. Each panel is one tissue type. The plot has `nSegments x linesPerSegment` phase encoding lines (default 3 x 64 = 192). The line direction alternates left-to-right then right-to-left, so it looks like a zigzag EPI readout. The brightness/color strength of each line comes directly from that tissue's simulated signal for the segment that acquired it — each block of `linesPerSegment` lines shows its own decay envelope restarting from a bright value at its segment's excitation, rather than one decayed train tiled identically across segments. Dotted horizontal guides mark the segment boundaries. This is useful because it shows how the WM and GM signals are actually placed into k-space, rather than just plotting signal versus echo number.

## Figure 6: 2D PSF k-space

<img width="1361" height="659" alt="Image" src="https://github.com/user-attachments/assets/96d6312e-dbfb-4a72-9cd9-d85d3bbb617a" />

Figure 6 shows the 2D PSF in k-space calculated from the filled 2D EPI trajectory. After the WM/GM signal train is weighted into k-space, the code examines the spatial frequency distribution to show how the echo train shapes the encoded signal before the inverse transform. This is useful because it reveals how strongly the EPI weighting tapers across the readout and whether the resulting k-space coverage preserves the expected signal distribution for each tissue.

## Figure 7: PSF

<img width="1177" height="786" alt="Image" src="https://github.com/user-attachments/assets/f11c65a5-cd4b-44cf-9fec-e91388b353b2" />

The Point Spread Function (PSF) represents how a point object maps into image space given the simulated k-space sampling and echo weighting from the EPI train. The PSF shows the main lobe (resolution) and sidelobes (ringing or aliasing) whose pattern depends on readout window width, partition sampling, and any apodisation due to echo amplitude envelope. Inspect the PSF for expected main lobe width consistent with k-space coverage and for sidelobe patterns that reflect sampling density and ordering. Unexpected asymmetry or artifact structures often indicate issues in k-space ordering or missing complex phase terms in k-space accumulation.

## Figure 8: TI Optimisation

<img width="1177" height="786" alt="Image" src="https://github.com/user-attachments/assets/ae690e5f-7415-4896-bebc-fefe192f4601" />

This figure presents the residual signal of the tissue of interest as `TI` varies. It is used to select the inversion time that minimises that tissue's signal. Key validation checks are: the TI sweep range covers the predicted null point, the step size is sufficiently fine to resolve the minimum, and the curve around the minimum is smooth (indicating stable numeric computation rather than noise). A well formed dip near the analytic null time indicates that inversion recovery and the `Mz_init` calculation are functioning as intended.

## Figure: Comparison

<img width="1446" height="938" alt="Image" src="https://github.com/user-attachments/assets/38756557-cea8-451e-8296-052a27f3f5ec" />

## Known limitations

- steady-state magnetization across shots. Every shot currently starts its FLAIR prep from full equilibrium (`Mz=1`), which is only strictly accurate for the first shot of a scan.

- segment-to-k-space ordering. Segments are currently assigned to *contiguous* blocks of ky lines. Real segmented EPI could lead to the interleaving of segments across k-space instead (segment 1 acquires every Nth line, spanning the full width of k-space) so that each shot alone gives a coarse full-FOV view. 

- Real scanners excite edge partitions slightly less accurately than centre partitions or any residual T1 aturation between shots.

- No partial Fourier acceleration option. These are the other major levers (besides segmentation) that real 3D-EPI FLAIR protocols use to keep echo trains short enough to leave usable signal by the end of k-space at high resolution.

- Off-resonance/B0-inhomogeneity-driven geometric distortion, a defining aspect of EPI, is not modelled here. `T2star` only affects amplitude decay here, not phase/spatial distortion.

`EPI_vs_TSE_comparison.m` simplifications: 
- The TSE readout uses a constant 180 degree CPMG refocusing train. Clinical 3D TSE or SPACE-FLAIR protocols use variable (ramped) refocusing flip angles well below 180 degrees, which this model does not represent.
- B1+ inhomogeneity —> one of the largest practical limitations at 7T for both readouts is not modelled at all.

---

## Getting Started with MATLAB

### How to Open and Use in MATLAB

1. **Clone or download the repository** to your local machine:
   ```bash
   git clone https://github.com/viilerr/3D-EPG-EPI-FLAIR.git
   cd 3D-EPG-EPI-FLAIR
   ```

2. **Open MATLAB** and navigate to the repository folder:
   ```matlab
   cd /path/to/3D-EPG-EPI-FLAIR
   ```

3. **Add the repository to your MATLAB path** (optional but recommended):
   ```matlab
   addpath(genpath(pwd))
   ```

4. **Run the test script** to visualize and validate the simulations:
   ```matlab
   Test_EPG_FLAIR_Explained
   ```
   This generates eight interactive figures showing FLAIR recovery, EPI readout, EPG state evolution, k-space trajectories, and point-spread functions.

5. **Use the core function** in your own code:
   ```matlab
   [F0,K3D,Zn,F,info] = EPG_FLAIR(deg2rad(90),0.8,T1,T2,3000, ...
       'sequence','epi3d','T2star',40,'T2prepTE',50, ...
       'nPE',64,'nPartitions',32);
   ```

### Requirements
- MATLAB R2018b or later
- No additional toolboxes required
