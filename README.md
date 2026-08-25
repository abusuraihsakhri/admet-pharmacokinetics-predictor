# ADMET & Pharmacokinetics Predictor

A high-performance, pure Python standard library computational engine for early-stage drug discovery, physicochemical property profiling, rule-based druglikeness filtering, Quantitative Estimate of Drug-likeness (QED), Central Nervous System Multiparameter Optimization (CNS MPO), Blood-Brain Barrier (BBB) permeability prediction, and analytical 1-/2-compartment pharmacokinetic (PK) simulation.

---

## Key Features

- **Druglikeness Rule Evaluation**:
  - **Lipinski's Rule of Five (Ro5)**: Molecular Weight $\le 500\text{ Da}$, $\text{LogP} \le 5.0$, $\text{HBD} \le 5$, $\text{HBA} \le 10$ (tolerates $\le 1$ violation).
  - **Veber's Bioavailability Filter**: Rotatable Bonds $\le 10$, $\text{TPSA} \le 140\text{ \AA}^2$.
  - **Egan's Egg Filter**: $-1.0 \le \text{LogP} \le 5.88$, $\text{TPSA} \le 131.6\text{ \AA}^2$.
  - **Ghose Filter**: $160 \le \text{MW} \le 480$, $-0.4 \le \text{LogP} \le 5.6$, $40 \le \text{MR} \le 130$, $20 \le \text{Atoms} \le 70$.
  - **Muegge (Bayer) Filter**: Structural and physicochemical criteria for drug-like chemical space.
  - **Lead-Likeness Filter**: $150 \le \text{MW} \le 350$, $\text{LogP} \le 3.0$, $\text{Rotatable Bonds} \le 7$.

- **Quantitative Scoring & Optimization**:
  - **QED (Bickerton et al., Nature Chem 2012)**: Asymmetric double sigmoidal desirability functions for 8 molecular descriptors aggregated via weighted geometric mean:
    $$\text{QED} = \exp\left(\frac{\sum_i w_i \ln d_i}{\sum_i w_i}\right)$$
  - **CNS MPO (Pfizer, Wager et al. 2010, 2016)**: Continuous multi-parameter optimization across 6 properties (ClogP, ClogD7.4, MW, TPSA, HBD, $\text{p}K_a$), predicting BBB permeability and brain exposure.
  - **LogBB Partition Model**: $\text{LogBB} = 0.152 \cdot \text{LogP} - 0.0148 \cdot \text{TPSA} + 0.139$.

- **In Vitro & In Vivo ADMET Predictions**:
  - **Absorption**: Human Intestinal Absorption ($\text{HIA } \%$), Caco-2 apparent permeability ($P_{\text{app}}$ in $10^{-6}\text{ cm/s}$), P-gp substrate & inhibitor risk.
  - **Distribution**: Plasma Protein Binding ($\text{PPB } \%$), Fraction unbound ($f_u$), Steady-State Volume of Distribution ($V_{d,ss}$ in $\text{L/kg}$).
  - **Metabolism**: Cytochrome P450 isoform inhibition liabilities (CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4).
  - **Excretion**: Total body clearance ($\text{mL/min/kg}$), elimination rate constant $k_e$, elimination half-life $t_{1/2}$.
  - **Toxicity Liabilities**: hERG cardiac potassium channel risk, Drug-Induced Liver Injury (DILI) alert, Ames mutagenicity screening.

- **Pharmacokinetic Simulation Engine**:
  - **Single Oral Dose**: Analytical Bateman function:
    $$C(t) = \frac{D \cdot F \cdot k_a}{V_d (k_a - k_e)} \left(e^{-k_e t} - e^{-k_a t}\right)$$
    $$T_{\max} = \frac{\ln(k_a / k_e)}{k_a - k_e}, \quad C_{\max} = C(T_{\max}), \quad \text{AUC}_{0-\infty} = \frac{D \cdot F}{\text{CL}}$$
  - **Intravenous Bolus**: Single-compartment exponential decay $C(t) = \frac{D}{V_d} e^{-k_e t}$.
  - **Multiple-Dose Oral (Steady State)**: Superposition model calculating accumulation ratio $R_{\text{acc}} = \frac{1}{1 - e^{-k_e \tau}}$ and average steady-state concentration $C_{ss,\text{avg}} = \frac{D \cdot F}{\text{CL} \cdot \tau}$.

---

## Installation & Setup

Zero external dependencies required; runs on pure standard library Python 3.9+.

```bash
git clone https://github.com/abusuraihsakhri/admet-pharmacokinetics-predictor.git
cd admet-pharmacokinetics-predictor
```

---

## CLI Usage

### 1. Interactive Studio
Launch the interactive wizard:
```bash
python cli.py
# or
python cli.py --interactive
```

### 2. Reference Benchmark Demonstration
Run built-in benchmark evaluations across reference therapeutics (Aspirin, Caffeine, Ibuprofen, Atorvastatin, Imatinib, Morphine, Vancomycin, Metformin, Diazepam):
```bash
python cli.py --demo
```

### 3. Evaluate Single Molecule Candidate
```bash
python cli.py evaluate --name "Lead-Candidate-42" \
  --mw 345.2 \
  --logp 2.4 \
  --hbd 1 \
  --hba 4 \
  --tpsa 58.5 \
  --rotatable 4 \
  --aromatic 2 \
  --heavy 24 \
  --mr 82.0 \
  --pka-base 8.4
```

### 4. Pharmacokinetic Simulation
Simulate oral concentration-time profiles:
```bash
python cli.py pk-sim --route oral --dose 100 --f 0.85 --ka 1.2 --ke 0.12 --vd 22.0 --duration 24.0
```

Simulate steady-state multiple dosing:
```bash
python cli.py pk-sim --route multi --dose 250 --f 0.80 --ka 1.1 --ke 0.08 --vd 30.0 --tau 12.0 --doses 7
```

### 5. Batch CSV Processing
Evaluate a collection of molecules in batch:
```bash
python cli.py batch --input molecules.csv --output results.csv
```

---

## Python API Example

```python
from admet_predictor import (
    MoleculeProperties,
    ADMETPredictor,
    PharmacokineticSimulator,
    REFERENCE_DRUGS,
)

# Define or look up a molecule
mol = MoleculeProperties(
    name="Imatinib Analog",
    mw=450.5,
    logp=3.1,
    hbd=2,
    hba=6,
    tpsa=78.2,
    rotatable_bonds=6,
    aromatic_rings=3,
    heavy_atoms=33,
    pka_base=8.0,
)

# Run full evaluation
report = ADMETPredictor.evaluate_candidate(mol)
print(f"Overall Score: {report.overall_druglikeness_score}/100")
print(f"QED Score: {report.qed.qed_score}")
print(f"CNS MPO: {report.cns_mpo.score} ({report.cns_mpo.cns_permeability_likelihood})")
print(f"Estimated Half-life: {report.admet_prediction.elimination_half_life_hr} hr")

# Simulate Pharmacokinetics
sim = PharmacokineticSimulator.simulate_oral_single(
    dose_mg=200.0,
    bioavailability_f=0.85,
    ka_hr=1.2,
    ke_hr=math.log(2) / report.admet_prediction.elimination_half_life_hr,
    vd_l=report.admet_prediction.vd_ss_l_kg * 70.0, # for 70 kg human
)
print(f"Peak Concentration Cmax: {sim.cmax_mg_l:.3f} mg/L at Tmax: {sim.tmax_hr:.1f} hr")
```

---

## Running Test Suite

Execute the 30-case unit test suite:

```bash
python -m unittest test_admet_predictor.py
```

---

## License

MIT License. Designed for medicinal chemistry, AI drug discovery, and translational pharmacology workflows.
