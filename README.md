# ADMET & Pharmacokinetics Predictor

### [Open the Live Application →](https://abusuraihsakhri.github.io/admet-pharmacokinetics-predictor/)

A Python tool for descriptor-based drug-likeness screening, heuristic ADMET estimates, CNS multiparameter optimization (CNS MPO) scoring, and one-compartment pharmacokinetic simulation.

The Python engine uses only the standard library. The browser application runs the same engine client-side with Pyodide.

## What it provides

- Rule-based filters: Lipinski, Veber, Egan, Ghose, Muegge, and lead-likeness.
- A simplified **QED-like** descriptor score.
- Six-parameter CNS MPO scoring using LogP, LogD7.4, molecular weight, TPSA, H-bond donors, and basic pKa.
- Descriptor-based heuristic estimates for HIA, Caco-2 permeability, P-gp, PPB, distribution volume, CYP inhibition flags, clearance, half-life, hERG, DILI, and Ames risk.
- One-compartment PK models for oral single dose, IV bolus, and repeated oral dosing.
- Analytical steady-state oral Cmax, Cmin, Cavg, Tmax, AUC per interval, and peak accumulation ratio.
- Batch CSV processing and reference-drug examples.
- Responsive browser UI with light/dark themes.

## Browser application

Open the live application above. The first load downloads the pinned Pyodide runtime from jsDelivr; calculations then execute in the browser using the repository's Python source.

The application has two workspaces:

1. **Molecule screen** — enter physicochemical descriptors and review rule-based filters, CNS MPO, the QED-like score, and heuristic ADMET flags.
2. **PK simulation** — simulate oral single-dose, IV bolus, or repeated oral dosing and inspect concentration-time output.

### Privacy

Form values are processed locally in the browser. The application does not submit molecule or PK inputs to this repository or to an application backend. Loading the page and Pyodide still makes normal network requests to GitHub Pages and the jsDelivr CDN.

## Important limitations

This repository is an exploratory screening and educational tool, not a validated ADMET prediction platform.

- The ADMET equations are transparent descriptor heuristics. They are not trained or externally validated predictive models.
- The reported **QED-like score is not canonical Bickerton QED**. Exact QED requires molecular-structure-derived inputs, including structural alerts and aromatic atom proportion, which this descriptor-only interface does not have.
- CNS MPO is a property-alignment score, not a direct measurement of blood-brain barrier permeability or clinical CNS exposure.
- If LogD7.4 is not supplied, the engine estimates it from LogP and one acidic or basic pKa; ampholytes and complex ionization behavior are not modeled rigorously.
- PK simulations use standard one-compartment analytical assumptions. Parameters such as F, ka, ke, and Vd must be supplied from appropriate evidence.
- Outputs should not be used for clinical, regulatory, dosing, or patient-care decisions.

## Command-line use

Requires Python 3.10 or newer.

```bash
git clone https://github.com/abusuraihsakhri/admet-pharmacokinetics-predictor.git
cd admet-pharmacokinetics-predictor
python -m pip install -e .
```

Evaluate a molecule:

```bash
admet-predictor evaluate \
  --name Candidate-01 \
  --mw 320.5 \
  --logp 2.8 \
  --hbd 2 \
  --hba 4 \
  --tpsa 60
```

JSON output:

```bash
admet-predictor evaluate --mw 320.5 --logp 2.8 --json
```

Reference-drug example:

```bash
admet-predictor ref Aspirin
```

Repeated oral dosing:

```bash
admet-predictor pk-sim \
  --route multi \
  --dose 250 \
  --f 0.8 \
  --ka 1.0 \
  --ke 0.1 \
  --vd 30 \
  --tau 12 \
  --doses 7
```

Batch CSV processing:

```bash
admet-predictor batch --input sample.csv --output results.csv
```

## Python API

```python
from admet_predictor import (
    ADMETPredictor,
    MoleculeProperties,
    PharmacokineticSimulator,
)

molecule = MoleculeProperties(
    name="Candidate-A",
    mw=325.4,
    logp=2.4,
    hbd=1,
    hba=4,
    tpsa=55.0,
    rotatable_bonds=4,
    aromatic_rings=2,
    heavy_atoms=23,
    molar_refractivity=78.0,
)

report = ADMETPredictor.evaluate_candidate(molecule)
print(report.overall_druglikeness_score)
print(report.cns_mpo.score)

simulation = PharmacokineticSimulator.simulate_oral_single(
    dose_mg=100,
    bioavailability_f=0.85,
    ka_hr=1.2,
    ke_hr=0.15,
    vd_l=25,
)
print(simulation.cmax_mg_l, simulation.tmax_hr)
```

## Local development

Install the package and test dependency:

```bash
python -m pip install -e . pytest
```

Run the test suite:

```bash
pytest -q
```

Run the browser application locally from the repository root:

```bash
python -m http.server 8000
```

Then open `http://localhost:8000/`. A local HTTP server is required because the browser app fetches the Python engine source.

## Technology

- Python 3.10+
- Standard-library Python runtime
- Pyodide 0.29.5 for in-browser Python execution
- HTML, CSS, and vanilla JavaScript
- GitHub Actions for CI and GitHub Pages deployment

GitHub Actions test Python 3.10, 3.11, 3.12, and 3.13. The web application targets current versions of Chrome, Edge, Firefox, and Safari with WebAssembly support.

## Method references

- Lipinski CA et al. *Adv Drug Deliv Rev.* 2001;46:3–26. doi: [10.1016/S0169-409X(00)00129-0](https://doi.org/10.1016/S0169-409X(00)00129-0)
- Veber DF et al. *J Med Chem.* 2002;45:2615–2623. doi: [10.1021/jm020017n](https://doi.org/10.1021/jm020017n)
- Wager TT et al. *ACS Chem Neurosci.* 2010;1:435–449. doi: [10.1021/cn100008c](https://doi.org/10.1021/cn100008c)
- Bickerton GR et al. *Nat Chem.* 2012;4:90–98. doi: [10.1038/nchem.1243](https://doi.org/10.1038/nchem.1243)

The references describe the underlying published concepts. They should not be interpreted as validation of this repository's heuristic ADMET equations or simplified QED-like implementation.

## License

MIT. See [LICENSE](LICENSE).
