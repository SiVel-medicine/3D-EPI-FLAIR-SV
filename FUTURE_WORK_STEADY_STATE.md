# Future Work: Steady-State Magnetization Chaining

## Status: not finished as I couldn't fully comprehend


Every shot starts its FLAIR preparation from full equilibrium magnetization (`opts.zinit`, default `Mz=1`). That's accurate for the *first* shot of a scan, but not for shot 2 onward: a real scanner doesn't let every tissue fully recover to quilibrium between 96 back-to-back shots. The true starting magnetization for most shots is somewhere between equilibrium and wherever the previous shot's readout left it. It is a genuine "steady state" that the current equilibrium-every-shot model doesn't capture.

TI is calibrated so CSF crosses zero magnetization starting from a specific pre-inversion Mz. If the real starting Mz differs shot-to-shot from what TI was calibrated for, CSF nulling quality will vary across the acquisition which seems to be an effect from a more realistic model.

## What was found (and why it's unresolved)

Chaining shots can be implemented. Carrying the magnetisation (Mz) from one shot to the next is easy to implement, and it reliably settles into a steady state. This was checked and gave me identical results in MATLAB, but it breaks the CSF nulling. At steady state, CSF's Mz sits far from full equilibrium, so the original TI of 3000 ms (tuned for starting at equilibrium) no longer nulls CSF.
I tried a fix which didn't actually work, by adding a 200-250 ms gap between shots which looked promising early on, but that was based on a rough estimate. Recalculated properly, the TI that would null CSF at true steady state is only about 250-500 ms, which is unrealistically short for FLAIR. This maybe due to the model being too simplified. It repeats the full FLAIR preparation (T2-prep + inversion + TI) before every shot. Real 3D sequences should share one inversion pulse across several partitions or segments before re-preparing, which is why they get good CSF nulling at realistic TIs and this model doesn't.

Tldr: The chaining code works correctly, but it basically revealed that the model's assumption of a full re-preparation before every shot may be slightly unrealistic. Fixing the code alone doesn't fix that assumption and needs to be more time and expertise to sort out.

## What real segmented IR sequences do about this

Papers such as Visser et al. (2010) describes a few practical 7T 3D-FLAIR implementations, with the general literature on segmented inversion-recovery sequences describing two standard mitigations:

- **Dummy prep shots** -> run a few throwaway shots (prep + excite + readout, discarded) before real data acquisition begins, so Mz reaches steady state before any real k-space lines are collected.
- **Ramped/variable flip angle**: gradually increase the excitation flip angle over the first several real shots instead of jumping straight to the target angle, smoothing the approach to steady state.

These do not by themselves fix a TI that's miscalibrated for the converged steady-state Mz.

## What code would be needed to implement this properly

If someone were to consider to pick this up:

- **Change the timing model first.** -> Let one inversion pulse serve many partitions or segments before the next T2-prep + inversion cycle. This also changes the shot-duration and total-scan-time calculations, not just the Mz chaining.

- **Recalibrate TI against the true steady-state Mz** 

- **Add a warm-up of dummy shots or a ramped file angle.**

- **Check things by confirming Mz converges and that csf is nulled properly.** 
