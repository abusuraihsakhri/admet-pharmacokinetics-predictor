# ADMET & Pharmacokinetics Predictor Engine

A Python computational chemistry, medicinal chemistry, and pharmacokinetics evaluation engine. Evaluates small molecule drug-likeness rules (Lipinski Rule of 5, Veber, Egan, Ghose, Muegge, Lead-likeness), quantitative drug-likeness (QED), Central Nervous System Multiparameter Optimization (CNS MPO), absorption/distribution/metabolism/excretion/toxicity (ADMET) risks, and one-compartment pharmacokinetic concentration-time simulations.

Requires Python standard library only (zero external runtime dependencies).

---

## Features

- **Drug-Likeness Rules & Filtering:**
  - **Lipinski Rule of 5:** MW $\le$ 500 Da, LogP $\le$ 5.0, HBD $\le$ 5, HBA $\le$ 10.
  - **Veber Filter:** Rotatable bonds $\le$ 10, TPSA $\le$ 140 $\text{\AA}^2$.
  - **Egan, Ghose, Muegge (Bayer), & Lead-Likeness Filters.**
- **Quantitative Estimation of Drug-Likeness (QED):** Bickerton et al. asymmetric desirability function across 8 physicochemical descriptors.
- **CNS Multiparameter Optimization (CNS MPO):** Wager et al. 6-parameter scoring function (0-6.0 scale) and predicted LogBB brain penetration likelihood.
- **In-Silico ADMET Property Profiling:**
  - **Absorption:** Human Intestinal Absorption (HIA %), Caco-2 permeability, P-gp substrate/inhibition risk.
  - **Distribution:** Plasma protein binding (PPB %), volume of distribution ($V_{d,\text{ss}}$).
  - **Metabolism:** CYP450 inhibition risk profiling (CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4).
  - **Excretion & Toxicity:** Intrinsic clearance, elimination half-life, hERG cardiotoxicity, DILI hepatotoxicity, Ames mutagenicity.
- **Pharmacokinetic (PK) Simulator:**
  - Single-dose oral absorption model (Bateman function).
  - IV bolus elimination model.
  - Multi-dose steady-state oral kinetics ($C_{\text{max}}$, $C_{\text{min}}$, $C_{\text{ss,avg}}$, accumulation ratio $R$).
- **Reference Drug Library:** Built-in benchmarking profiles for Aspirin, Caffeine, Ibuprofen, Atorvastatin, Imatinib, Morphine, Vancomycin, Metformin, and Diazepam.
- **Batch CSV Processing:** High-throughput candidate library evaluation.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies. `pytest` is optional for running tests.

```bash
git clone https://github.com/abusuraihsakhri/admet-pharmacokinetics-predictor.git
cd admet-pharmacokinetics-predictor
```

---

## CLI Usage

### 1. Evaluate Single Molecule
```bash
python cli.py evaluate --name Candidate-01 --mw 320.5 --logp 2.8 --hbd 2 --hba 4 --tpsa 60.0
```
Output as JSON:
```bash
python cli.py evaluate --name Candidate-01 --mw 320.5 --logp 2.8 --json
```

### 2. Reference Drug Benchmark
```bash
python cli.py ref Aspirin
python cli.py ref Imatinib --json
python cli.py --demo
```

### 3. Pharmacokinetic Simulation
Simulate oral single dose:
```bash
python cli.py pk-sim --route oral --dose 100 --f 0.85 --ka 1.2 --ke 0.15 --vd 25
```
Simulate multi-dose regimen:
```bash
python cli.py pk-sim --route multi --dose 250 --f 0.80 --ka 1.0 --ke 0.10 --vd 30 --tau 12 --doses 7 --json
```

### 4. Batch CSV Evaluation
```bash
python cli.py batch --input sample.csv --output results.csv
```

---

## Python API Quickstart

```python
from admet_predictor import MoleculeProperties, ADMETPredictor, PharmacokineticSimulator

# 1. Evaluate small molecule properties
mol = MoleculeProperties(
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

report = ADMETPredictor.evaluate_candidate(mol)
print(f"Lipinski Pass: {report.lipinski.passes} (Violations: {report.lipinski.violations})")
print(f"QED Score: {report.qed.qed_score:.3f} [{report.qed.druglikeness_grade}]")
print(f"CNS MPO Score: {report.cns_mpo.score:.2f} [{report.cns_mpo.cns_permeability_likelihood}]")
print(f"Overall Score: {report.overall_druglikeness_score:.1f} / 100")

# 2. Simulate single-dose oral PK profile
sim = PharmacokineticSimulator.simulate_oral_single(
    dose_mg=100.0,
    bioavailability_f=0.85,
    ka_hr=1.2,
    ke_hr=0.15,
    vd_l=25.0,
)
print(f"Cmax: {sim.cmax_mg_l:.3f} mg/L at Tmax: {sim.tmax_hr:.2f} hr | t1/2: {sim.half_life_hr:.2f} hr")
```

---

## Running Tests

Run the test suite using standard `unittest` or `pytest`:

```bash
python test_admet_predictor.py
# or
pytest -v
```

