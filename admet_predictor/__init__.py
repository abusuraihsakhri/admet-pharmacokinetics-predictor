#!/usr/bin/env python3
"""
ADMET & Pharmacokinetics Predictor
===================================
A comprehensive, pure standard library engine for physicochemical property evaluation,
drug-likeness filtering (Lipinski Rule of 5, Veber, Egan, Ghose, Muegge, Lead-likeness),
Quantitative Estimate of Drug-likeness (QED, Bickerton et al. 2012),
CNS Multiparameter Optimization (CNS MPO, Wager et al. 2016),
Blood-Brain Barrier (BBB) Permeability, ADMET risk profiling,
and 1-compartment / 2-compartment pharmacokinetic (PK) simulation.
"""

from __future__ import annotations

import math
import json
import csv
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any, Union

__version__ = "2.0.0"

# ============================================================================
# Data Models
# ============================================================================

@dataclass
class MoleculeProperties:
    """Physicochemical descriptors of a drug candidate."""
    name: str = "Candidate"
    mw: float = 300.0  # Molecular weight (g/mol or Da)
    logp: float = 2.0  # Octanol-water partition coefficient (ClogP / ALogP)
    hbd: int = 2       # Hydrogen bond donors (OH, NH count)
    hba: int = 4       # Hydrogen bond acceptors (N, O count)
    tpsa: float = 60.0 # Topological Polar Surface Area (Å²)
    rotatable_bonds: int = 4 # Number of freely rotatable bonds
    aromatic_rings: int = 2  # Number of aromatic rings
    heavy_atoms: int = 22    # Number of non-hydrogen atoms
    molar_refractivity: float = 75.0 # Molar refractivity (MR)
    pka_base: Optional[float] = None # Most basic pKa (if applicable)
    pka_acid: Optional[float] = None # Most acidic pKa (if applicable)
    logd74: Optional[float] = None   # Distribution coefficient at physiological pH 7.4
    fsp3: float = 0.35 # Fraction of sp3 hybridized carbons (0-1)
    charge: int = 0    # Net formal charge at pH 7.4

    def __post_init__(self):
        if self.mw <= 0:
            raise ValueError(f"Molecular weight must be positive, got {self.mw}")
        if self.hbd < 0 or self.hba < 0 or self.rotatable_bonds < 0:
            raise ValueError("HBD, HBA, and rotatable bonds count cannot be negative")
        if self.tpsa < 0:
            raise ValueError("TPSA cannot be negative")
        if self.heavy_atoms < 0:
            raise ValueError("Heavy atom count cannot be negative")
        if not (0.0 <= self.fsp3 <= 1.0):
            self.fsp3 = max(0.0, min(1.0, self.fsp3))
        # Default LogD7.4 approximation if not provided
        if self.logd74 is None:
            if self.pka_base is not None and self.pka_base > 7.4:
                ionization_ratio = 1.0 + 10 ** (self.pka_base - 7.4)
                self.logd74 = round(self.logp - math.log10(ionization_ratio), 2)
            elif self.pka_acid is not None and self.pka_acid < 7.4:
                ionization_ratio = 1.0 + 10 ** (7.4 - self.pka_acid)
                self.logd74 = round(self.logp - math.log10(ionization_ratio), 2)
            else:
                self.logd74 = self.logp


@dataclass
class FilterResult:
    """Result of a specific drug-likeness filter."""
    rule_name: str
    passes: bool
    violations: int
    criteria_results: Dict[str, bool]
    rationale: str
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class QEDResult:
    """Quantitative Estimate of Drug-likeness (Bickerton et al. 2012)."""
    qed_score: float  # Scale 0.0 to 1.0
    druglikeness_grade: str # High (>0.67), Moderate (0.49-0.67), Low (<0.49)
    individual_desirabilities: Dict[str, float]
    weights: Dict[str, float]


@dataclass
class CNSMPOResult:
    """Pfizer CNS Multiparameter Optimization (Wager et al. 2010, 2016)."""
    score: float # Scale 0.0 to 6.0
    cns_permeability_likelihood: str # High (>=4.0), Moderate (2.5-4.0), Low (<2.5)
    logbb_pred: float # Estimated Log([Brain]/[Plasma])
    component_scores: Dict[str, float]
    risk_factors: List[str]


@dataclass
class ADMETPropertyEstimate:
    """Individual ADMET parameter predictions."""
    hia_pct: float            # Human Intestinal Absorption (%)
    caco2_perm_cm_s: float    # Caco-2 permeability (10^-6 cm/s)
    caco2_class: str          # High, Moderate, Low
    pgp_substrate_risk: str   # High, Moderate, Low
    pgp_inhibitor_risk: str   # High, Moderate, Low
    ppb_pct: float            # Plasma Protein Binding (%)
    fraction_unbound: float   # fu (0 to 1)
    vd_ss_l_kg: float         # Volume of distribution steady-state (L/kg)
    vd_classification: str    # Low (<0.7 L/kg), Moderate (0.7-2 L/kg), High (>2 L/kg)
    cyp_inhibitions: Dict[str, str] # CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4 -> Low/Medium/High
    clearance_ml_min_kg: float # Total clearance estimate (mL/min/kg)
    elimination_half_life_hr: float # Half-life t1/2 (hours)
    herg_risk: str            # High, Moderate, Low (cardiotoxicity)
    dili_risk: str            # High, Moderate, Low (hepatotoxicity)
    ames_mutagenicity: str    # Positive / Negative likelihood


@dataclass
class PKSimulationPoint:
    """Time-series point for PK simulation."""
    time_hr: float
    plasma_conc_mg_l: float
    central_conc_mg_l: Optional[float] = None
    peripheral_conc_mg_l: Optional[float] = None


@dataclass
class PKSimulationResult:
    """Comprehensive pharmacokinetic simulation result."""
    dosing_route: str # IV_BOLUS, ORAL_SINGLE, ORAL_MULTIPLE
    dose_mg: float
    bioavailability_f: float
    half_life_hr: float
    elimination_rate_ke: float
    absorption_rate_ka: Optional[float]
    volume_distribution_l: float
    clearance_l_hr: float
    cmax_mg_l: float
    tmax_hr: float
    auc_0_inf_mg_hr_l: float
    auc_tau_mg_hr_l: Optional[float] = None
    c_ss_avg_mg_l: Optional[float] = None
    c_ss_min_mg_l: Optional[float] = None
    c_ss_max_mg_l: Optional[float] = None
    accumulation_ratio: Optional[float] = None
    concentration_curve: List[PKSimulationPoint] = field(default_factory=list)


@dataclass
class ComprehensiveADMETReport:
    """Full comprehensive evaluation report for a drug candidate."""
    molecule: MoleculeProperties
    lipinski: FilterResult
    veber: FilterResult
    egan: FilterResult
    ghose: FilterResult
    muegge: FilterResult
    lead_likeness: FilterResult
    qed: QEDResult
    cns_mpo: CNSMPOResult
    admet_prediction: ADMETPropertyEstimate
    overall_druglikeness_score: float # 0 to 100
    overall_assessment: str
    recommendations: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


# ============================================================================
# Rule Evaluators
# ============================================================================

class LipinskiRuleOf5:
    """Lipinski's Rule of 5: MW<=500, LogP<=5, HBD<=5, HBA<=10."""
    MW_MAX = 500.0
    LOGP_MAX = 5.0
    HBD_MAX = 5
    HBA_MAX = 10

    @classmethod
    def evaluate(cls, mol: MoleculeProperties) -> FilterResult:
        criteria = {
            "MW <= 500 Da": mol.mw <= cls.MW_MAX,
            "LogP <= 5.0": mol.logp <= cls.LOGP_MAX,
            "HBD <= 5": mol.hbd <= cls.HBD_MAX,
            "HBA <= 10": mol.hba <= cls.HBA_MAX,
        }
        violations = sum(1 for passed in criteria.values() if not passed)
        passes = violations <= 1

        failed = [k for k, v in criteria.items() if not v]
        if violations == 0:
            rationale = "Complies fully with all 4 Lipinski criteria (0 violations)."
        elif violations == 1:
            rationale = f"Complies with Lipinski Rule of 5 (1 permitted violation: {failed[0]})."
        else:
            rationale = f"Fails Lipinski Rule of 5 with {violations} violations: {', '.join(failed)}."

        return FilterResult(
            rule_name="Lipinski Rule of 5",
            passes=passes,
            violations=violations,
            criteria_results=criteria,
            rationale=rationale,
            details={"max_allowed_violations": 1}
        )


class VeberRule:
    """Veber et al. (2002): Rotatable bonds <= 10 and TPSA <= 140 Å²."""
    ROTATABLE_MAX = 10
    TPSA_MAX = 140.0

    @classmethod
    def evaluate(cls, mol: MoleculeProperties) -> FilterResult:
        criteria = {
            "Rotatable Bonds <= 10": mol.rotatable_bonds <= cls.ROTATABLE_MAX,
            "TPSA <= 140 Å²": mol.tpsa <= cls.TPSA_MAX,
        }
        violations = sum(1 for passed in criteria.values() if not passed)
        passes = violations == 0

        failed = [k for k, v in criteria.items() if not v]
        if passes:
            rationale = "Complies with Veber oral bioavailability criteria."
        else:
            rationale = f"Fails Veber criteria ({violations} violation(s): {', '.join(failed)})."

        return FilterResult(
            rule_name="Veber Filter",
            passes=passes,
            violations=violations,
            criteria_results=criteria,
            rationale=rationale
        )


class EganRule:
    """Egan et al. (2000) Pharmaco-kinetic filter: -1.0 <= LogP <= 5.88 and TPSA <= 131.6 Å²."""
    LOGP_MIN = -1.0
    LOGP_MAX = 5.88
    TPSA_MAX = 131.6

    @classmethod
    def evaluate(cls, mol: MoleculeProperties) -> FilterResult:
        criteria = {
            "-1.0 <= LogP <= 5.88": (cls.LOGP_MIN <= mol.logp <= cls.LOGP_MAX),
            "TPSA <= 131.6 Å²": mol.tpsa <= cls.TPSA_MAX,
        }
        violations = sum(1 for passed in criteria.values() if not passed)
        passes = violations == 0

        failed = [k for k, v in criteria.items() if not v]
        rationale = "Passes Egan egg filter for absorption" if passes else f"Fails Egan filter: {', '.join(failed)}"

        return FilterResult(
            rule_name="Egan Filter",
            passes=passes,
            violations=violations,
            criteria_results=criteria,
            rationale=rationale
        )


class GhoseFilter:
    """Ghose et al. (1999): 160<=MW<=480, -0.4<=LogP<=5.6, 40<=MR<=130, 20<=Atoms<=70."""
    MW_MIN, MW_MAX = 160.0, 480.0
    LOGP_MIN, LOGP_MAX = -0.4, 5.6
    MR_MIN, MR_MAX = 40.0, 130.0
    ATOMS_MIN, ATOMS_MAX = 20, 70

    @classmethod
    def evaluate(cls, mol: MoleculeProperties) -> FilterResult:
        criteria = {
            "160 <= MW <= 480 Da": cls.MW_MIN <= mol.mw <= cls.MW_MAX,
            "-0.4 <= LogP <= 5.6": cls.LOGP_MIN <= mol.logp <= cls.LOGP_MAX,
            "40 <= MR <= 130": cls.MR_MIN <= mol.molar_refractivity <= cls.MR_MAX,
            "20 <= Heavy Atoms <= 70": cls.ATOMS_MIN <= mol.heavy_atoms <= cls.ATOMS_MAX,
        }
        violations = sum(1 for passed in criteria.values() if not passed)
        passes = violations == 0
        failed = [k for k, v in criteria.items() if not v]
        rationale = "Passes Ghose filter" if passes else f"Fails Ghose filter ({violations} violations: {', '.join(failed)})"

        return FilterResult(
            rule_name="Ghose Filter",
            passes=passes,
            violations=violations,
            criteria_results=criteria,
            rationale=rationale
        )


class MueggeFilter:
    """Muegge / Bayer Filter (2001): 200<=MW<=600, -2<=LogP<=5, TPSA<=150, Rings<=7, Carbon>4, Hetero>1, Rot<=15, HBD<=5, HBA<=10."""
    @classmethod
    def evaluate(cls, mol: MoleculeProperties) -> FilterResult:
        criteria = {
            "200 <= MW <= 600 Da": 200.0 <= mol.mw <= 600.0,
            "-2.0 <= LogP <= 5.0": -2.0 <= mol.logp <= 5.0,
            "TPSA <= 150.0 Å²": mol.tpsa <= 150.0,
            "Aromatic/Rings <= 7": mol.aromatic_rings <= 7,
            "Rotatable Bonds <= 15": mol.rotatable_bonds <= 15,
            "HBD <= 5": mol.hbd <= 5,
            "HBA <= 10": mol.hba <= 10,
        }
        violations = sum(1 for passed in criteria.values() if not passed)
        passes = violations == 0
        failed = [k for k, v in criteria.items() if not v]
        rationale = "Passes Muegge (Bayer) filter" if passes else f"Fails Muegge filter: {', '.join(failed)}"

        return FilterResult(
            rule_name="Muegge Filter",
            passes=passes,
            violations=violations,
            criteria_results=criteria,
            rationale=rationale
        )


class LeadLikenessFilter:
    """Teague et al. (1999) Lead-likeness: MW 150-350 Da, LogP <= 3.0, Rotatable Bonds <= 7."""
    @classmethod
    def evaluate(cls, mol: MoleculeProperties) -> FilterResult:
        criteria = {
            "150 <= MW <= 350 Da": 150.0 <= mol.mw <= 350.0,
            "LogP <= 3.0": mol.logp <= 3.0,
            "Rotatable Bonds <= 7": mol.rotatable_bonds <= 7,
        }
        violations = sum(1 for passed in criteria.values() if not passed)
        passes = violations == 0
        failed = [k for k, v in criteria.items() if not v]
        rationale = "Candidate possesses ideal lead-like properties" if passes else f"Non-lead-like: {', '.join(failed)}"

        return FilterResult(
            rule_name="Lead-likeness Filter",
            passes=passes,
            violations=violations,
            criteria_results=criteria,
            rationale=rationale
        )


# ============================================================================
# Quantitative Estimate of Drug-likeness (QED)
# ============================================================================

class QEDCalculator:
    """
    Quantitative Estimate of Drug-likeness (Bickerton et al., Nature Chemistry 2012).
    Calculates asymmetric double sigmoidal desirability functions for 8 molecular descriptors:
    d_i(x) = a + b / (1 + exp(-(x - c + d/2)/e)) * (1 - 1 / (1 + exp(-(x - c - d/2)/f)))
    QED = exp( (sum_i w_i * ln d_i) / (sum_i w_i) )
    """
    ADS_PARAMS = {
        "mw": (0.47, 0.53, 345.96, 218.42, 47.93, 76.54, 0.99),
        "logp": (0.33, 0.67, 2.77, 3.49, 0.94, 1.25, 0.98),
        "hbd": (0.15, 0.85, 0.80, 2.76, 0.50, 0.80, 0.97),
        "hba": (0.20, 0.80, 3.60, 4.30, 0.80, 1.20, 0.98),
        "tpsa": (0.24, 0.76, 75.00, 78.00, 18.00, 26.00, 0.98),
        "rotatable": (0.20, 0.80, 4.00, 5.50, 1.20, 2.00, 0.98),
        "aromatic": (0.10, 0.90, 2.00, 2.20, 0.60, 1.00, 0.96),
        "alerts": (0.00, 1.00, 0.00, 0.50, 0.10, 0.50, 1.00),
    }

    WEIGHTS = {
        "mw": 0.66,
        "logp": 0.46,
        "hbd": 0.61,
        "hba": 0.05,
        "tpsa": 0.06,
        "rotatable": 0.65,
        "aromatic": 0.48,
        "alerts": 0.64,
    }

    @classmethod
    def _ads_fun(cls, x: float, a: float, b: float, c: float, d: float, e: float, f: float, dmax: float) -> float:
        try:
            exp1 = math.exp(-(x - c + d / 2.0) / e) if e != 0 else 0
            exp2 = math.exp(-(x - c - d / 2.0) / f) if f != 0 else 0
            s1 = 1.0 / (1.0 + exp1)
            s2 = 1.0 - 1.0 / (1.0 + exp2)
            val = (a + b * s1 * s2) / dmax
            return max(0.001, min(1.0, val))
        except OverflowError:
            return 0.001

    @classmethod
    def calculate(cls, mol: MoleculeProperties, alert_count: int = 0) -> QEDResult:
        desirabilities = {}
        raw_vals = {
            "mw": mol.mw,
            "logp": mol.logp,
            "hbd": float(mol.hbd),
            "hba": float(mol.hba),
            "tpsa": mol.tpsa,
            "rotatable": float(mol.rotatable_bonds),
            "aromatic": float(mol.aromatic_rings),
            "alerts": float(alert_count),
        }

        for prop, (a, b, c, d, e, f, dmax) in cls.ADS_PARAMS.items():
            val = raw_vals[prop]
            desirabilities[prop] = round(cls._ads_fun(val, a, b, c, d, e, f, dmax), 4)

        total_w = sum(cls.WEIGHTS.values())
        sum_weighted_log = sum(cls.WEIGHTS[k] * math.log(desirabilities[k]) for k in cls.WEIGHTS)
        qed = math.exp(sum_weighted_log / total_w)
        qed_score = round(max(0.0, min(1.0, qed)), 3)

        if qed_score >= 0.67:
            grade = "High Drug-Likeness"
        elif qed_score >= 0.49:
            grade = "Moderate Drug-Likeness"
        else:
            grade = "Low Drug-Likeness"

        return QEDResult(
            qed_score=qed_score,
            druglikeness_grade=grade,
            individual_desirabilities=desirabilities,
            weights=cls.WEIGHTS,
        )


# ============================================================================
# CNS MPO & Blood-Brain Barrier (BBB) Permeability
# ============================================================================

class CNSMPOPredictor:
    """
    Pfizer CNS Multiparameter Optimization (CNS MPO) scoring algorithm (Wager et al. 2010, 2016).
    Evaluates 6 properties with continuous monotonic/triangular ramps in [0, 1].
    """
    @classmethod
    def _ramp(cls, x: float, low: float, high: float, increasing: bool = False) -> float:
        if increasing:
            if x <= low: return 0.0
            if x >= high: return 1.0
            return (x - low) / (high - low)
        else:
            if x <= low: return 1.0
            if x >= high: return 0.0
            return (high - x) / (high - low)

    @classmethod
    def _tpsa_score(cls, tpsa: float) -> float:
        if 40.0 <= tpsa <= 90.0:
            return 1.0
        elif tpsa < 40.0:
            if tpsa <= 20.0: return 0.0
            return (tpsa - 20.0) / 20.0
        else:
            if tpsa >= 120.0: return 0.0
            return (120.0 - tpsa) / 30.0

    @classmethod
    def calculate(cls, mol: MoleculeProperties) -> CNSMPOResult:
        pka_val = mol.pka_base if mol.pka_base is not None else 7.0
        logd_val = mol.logd74 if mol.logd74 is not None else mol.logp

        s_clogp = cls._ramp(mol.logp, 3.0, 5.0, increasing=False)
        s_clogd = cls._ramp(logd_val, 2.0, 4.0, increasing=False)
        s_mw = cls._ramp(mol.mw, 360.0, 500.0, increasing=False)
        s_tpsa = cls._tpsa_score(mol.tpsa)
        s_hbd = cls._ramp(float(mol.hbd), 0.5, 3.5, increasing=False)
        s_pka = cls._ramp(pka_val, 8.0, 10.0, increasing=False)

        scores = {
            "ClogP": round(s_clogp, 3),
            "ClogD7.4": round(s_clogd, 3),
            "MW": round(s_mw, 3),
            "TPSA": round(s_tpsa, 3),
            "HBD": round(s_hbd, 3),
            "pKa_basic": round(s_pka, 3),
        }
        total_mpo = round(sum(scores.values()), 2)
        logbb = round(0.152 * mol.logp - 0.0148 * mol.tpsa + 0.139, 3)

        risk_factors = []
        if mol.tpsa > 90.0:
            risk_factors.append(f"High TPSA ({mol.tpsa:.1f} Å² > 90 Å²) reduces BBB penetration.")
        if mol.mw > 400.0:
            risk_factors.append(f"Molecular weight ({mol.mw:.1f} Da > 400 Da) restricts passive membrane diffusion.")
        if mol.hbd > 2:
            risk_factors.append(f"Excessive hydrogen bond donors (HBD={mol.hbd} > 2) impair brain uptake.")
        if mol.logp > 4.5:
            risk_factors.append(f"High lipophilicity (LogP={mol.logp:.2f}) increases nonspecific tissue binding.")
        if pka_val > 9.5:
            risk_factors.append(f"High basicity (pKa={pka_val:.1f}) increases lysosomal trapping.")

        if total_mpo >= 4.0 and logbb >= -0.5:
            likelihood = "High (CNS+)"
        elif total_mpo >= 2.5 and logbb >= -1.0:
            likelihood = "Moderate (CNS+/-)"
        else:
            likelihood = "Low / Non-permeant (CNS-)"

        return CNSMPOResult(
            score=total_mpo,
            cns_permeability_likelihood=likelihood,
            logbb_pred=logbb,
            component_scores=scores,
            risk_factors=risk_factors,
        )


# ============================================================================
# Comprehensive ADMET Predictor
# ============================================================================

class ADMETPredictor:
    """Predicts in vitro and in vivo ADMET parameters based on molecular descriptors."""

    @classmethod
    def predict_properties(cls, mol: MoleculeProperties) -> ADMETPropertyEstimate:
        # HIA %
        try:
            hia_exponent = 0.015 * mol.tpsa - 0.2 * mol.logp - 0.5
            hia_pct = round(100.0 / (1.0 + 10.0 ** hia_exponent), 1)
            hia_pct = max(5.0, min(99.0, hia_pct))
        except OverflowError:
            hia_pct = 10.0

        # Caco-2
        log_papp = -4.36 - 0.0104 * mol.tpsa + 0.32 * min(5.0, mol.logp) - 0.002 * min(600.0, mol.mw)
        caco2_perm = round(10.0 ** (log_papp + 6.0), 2)
        if caco2_perm >= 10.0:
            caco2_class = "High Permeability (>10x10^-6 cm/s)"
        elif caco2_perm >= 2.0:
            caco2_class = "Moderate Permeability (2-10x10^-6 cm/s)"
        else:
            caco2_class = "Low Permeability (<2x10^-6 cm/s)"

        # P-gp
        pgp_sub = "High (Likely P-gp substrate)" if (mol.mw > 400 and mol.hbd + mol.hba > 8) else ("Moderate" if mol.mw > 350 else "Low")
        pgp_inh = "High (Potential P-gp inhibitor)" if (mol.logp > 3.0 and mol.aromatic_rings >= 2) else ("Moderate" if mol.logp > 2.0 else "Low")

        # PPB & fu
        ppb_val = 100.0 / (1.0 + 10.0 ** (-0.7 * (mol.logp - 0.5) + 0.004 * mol.tpsa))
        ppb_pct = round(max(10.0, min(99.9, ppb_val)), 1)
        fu = round(max(0.001, min(1.0, (100.0 - ppb_pct) / 100.0)), 4)

        # Vd_ss
        vd_base = 0.5 + 0.8 * math.exp(max(-2.0, min(3.0, mol.logp * 0.4))) * (fu ** 0.5)
        if mol.pka_base is not None and mol.pka_base > 7.4:
            vd_base *= (1.0 + 0.3 * (mol.pka_base - 7.4))
        vd_ss = round(max(0.1, min(25.0, vd_base)), 2)
        if vd_ss < 0.7:
            vd_class = "Low (<0.7 L/kg, extracellular confinement)"
        elif vd_ss <= 2.0:
            vd_class = "Moderate (0.7-2.0 L/kg, total body water distribution)"
        else:
            vd_class = "High (>2.0 L/kg, extensive tissue distribution)"

        # CYP450
        cyps = {
            "CYP1A2": "High" if (mol.aromatic_rings >= 2 and mol.tpsa < 60 and mol.logp > 1.5) else ("Moderate" if mol.aromatic_rings >= 1 else "Low"),
            "CYP2C9": "High" if (mol.pka_acid is not None and mol.pka_acid < 6.0 and mol.logp > 2.0) else ("Moderate" if mol.logp > 2.5 else "Low"),
            "CYP2C19": "Moderate" if (mol.logp > 2.0 and mol.hba >= 3) else "Low",
            "CYP2D6": "High" if (mol.pka_base is not None and mol.pka_base >= 8.0 and mol.aromatic_rings >= 1) else ("Moderate" if mol.hbd >= 1 else "Low"),
            "CYP3A4": "High" if (mol.mw > 400 and mol.logp > 3.0 and mol.rotatable_bonds >= 4) else ("Moderate" if (mol.mw > 300 or mol.logp > 2.0) else "Low"),
        }

        # CL & t1/2
        cl_hep = max(0.5, min(20.0, 2.5 * (10 ** (0.25 * max(-1.0, min(3.0, mol.logp)))) * fu))
        cl_ren = max(0.1, min(10.0, 4.0 / (1.0 + 10 ** (0.4 * mol.logp)))) if mol.mw < 400 else 0.5
        cl_total = round(cl_hep + cl_ren, 2)
        cl_l_hr_kg = cl_total * 0.06
        ke = cl_l_hr_kg / vd_ss
        half_life = round(math.log(2.0) / max(0.001, ke), 2)
        half_life = max(0.2, min(72.0, half_life))

        # Toxicity liabilities
        herg_score = (1 if mol.logp > 3.0 else 0) + (1 if mol.aromatic_rings >= 2 else 0) + (2 if mol.pka_base is not None and mol.pka_base > 7.4 else 0)
        herg_risk = "High" if herg_score >= 3 else ("Moderate" if herg_score >= 2 else "Low")

        dili_score = (1 if mol.logp > 3.0 else 0) + (1 if mol.aromatic_rings >= 3 else 0) + (1 if mol.mw > 450 else 0)
        dili_risk = "High" if dili_score >= 3 else ("Moderate" if dili_score >= 1 else "Low")

        ames = "Positive Risk" if (mol.aromatic_rings >= 3 and mol.hbd >= 3 and mol.mw > 400) else "Negative (Low Risk)"

        return ADMETPropertyEstimate(
            hia_pct=hia_pct,
            caco2_perm_cm_s=caco2_perm,
            caco2_class=caco2_class,
            pgp_substrate_risk=pgp_sub,
            pgp_inhibitor_risk=pgp_inh,
            ppb_pct=ppb_pct,
            fraction_unbound=fu,
            vd_ss_l_kg=vd_ss,
            vd_classification=vd_class,
            cyp_inhibitions=cyps,
            clearance_ml_min_kg=cl_total,
            elimination_half_life_hr=half_life,
            herg_risk=herg_risk,
            dili_risk=dili_risk,
            ames_mutagenicity=ames,
        )

    @classmethod
    def evaluate_candidate(cls, mol: MoleculeProperties) -> ComprehensiveADMETReport:
        lipinski = LipinskiRuleOf5.evaluate(mol)
        veber = VeberRule.evaluate(mol)
        egan = EganRule.evaluate(mol)
        ghose = GhoseFilter.evaluate(mol)
        muegge = MueggeFilter.evaluate(mol)
        lead = LeadLikenessFilter.evaluate(mol)
        qed = QEDCalculator.calculate(mol)
        cns = CNSMPOPredictor.calculate(mol)
        admet = cls.predict_properties(mol)

        rule_score = (
            (1.0 if lipinski.passes else 0.5) * 30.0 +
            (1.0 if veber.passes else 0.0) * 20.0 +
            (1.0 if egan.passes else 0.0) * 10.0 +
            (1.0 if ghose.passes else 0.0) * 10.0 +
            (1.0 if muegge.passes else 0.0) * 10.0 +
            qed.qed_score * 20.0
        )
        overall_score = round(max(0.0, min(100.0, rule_score)), 1)

        recommendations = []
        if not lipinski.passes:
            recommendations.append(f"Address Lipinski violations ({lipinski.violations} found) to optimize oral absorption.")
        if not veber.passes:
            recommendations.append("Reduce rotatable bonds or polar surface area (TPSA) to improve permeability.")
        if admet.herg_risk == "High":
            recommendations.append("High hERG liability detected: consider decreasing basicity (pKa) or lipophilicity (LogP).")
        if admet.dili_risk == "High":
            recommendations.append("Elevated DILI hepatotoxicity alert: evaluate structural alerts and reactive metabolites.")
        if admet.cyp_inhibitions.get("CYP3A4") == "High":
            recommendations.append("Strong CYP3A4 inhibition predicted: investigate potential clinical drug-drug interactions (DDI).")
        if overall_score >= 80.0:
            assessment = "Excellent drug-like candidate with favorable physicochemical and ADMET profile."
        elif overall_score >= 60.0:
            assessment = "Acceptable drug candidate; target specific ADMET/physicochemical optimization."
        else:
            assessment = "Sub-optimal drug candidate with significant developability or ADMET liabilities."

        return ComprehensiveADMETReport(
            molecule=mol,
            lipinski=lipinski,
            veber=veber,
            egan=egan,
            ghose=ghose,
            muegge=muegge,
            lead_likeness=lead,
            qed=qed,
            cns_mpo=cns,
            admet_prediction=admet,
            overall_druglikeness_score=overall_score,
            overall_assessment=assessment,
            recommendations=recommendations,
        )


# ============================================================================
# Pharmacokinetic Simulation Engine
# ============================================================================

class PharmacokineticSimulator:
    """1-compartment and 2-compartment pharmacokinetic simulator."""

    @classmethod
    def simulate_oral_single(
        cls,
        dose_mg: float,
        bioavailability_f: float,
        ka_hr: float,
        ke_hr: float,
        vd_l: float,
        duration_hr: float = 24.0,
        num_points: int = 100,
    ) -> PKSimulationResult:
        if dose_mg <= 0 or bioavailability_f <= 0 or ka_hr <= 0 or ke_hr <= 0 or vd_l <= 0:
            raise ValueError("All PK parameters (dose, F, ka, ke, Vd) must be strictly positive.")

        if abs(ka_hr - ke_hr) < 1e-6:
            ka_hr += 1e-4

        cl_l_hr = ke_hr * vd_l
        t_half = math.log(2.0) / ke_hr
        t_max = math.log(ka_hr / ke_hr) / (ka_hr - ke_hr)
        t_max = max(0.0, t_max)

        def conc_func(t: float) -> float:
            c = (dose_mg * bioavailability_f * ka_hr / (vd_l * (ka_hr - ke_hr))) * (
                math.exp(-ke_hr * t) - math.exp(-ka_hr * t)
            )
            return max(0.0, c)

        c_max = conc_func(t_max)
        auc_inf = (dose_mg * bioavailability_f) / cl_l_hr

        dt = duration_hr / max(1, num_points - 1)
        curve = []
        for i in range(num_points):
            t = i * dt
            curve.append(PKSimulationPoint(time_hr=round(t, 2), plasma_conc_mg_l=round(conc_func(t), 4)))

        return PKSimulationResult(
            dosing_route="ORAL_SINGLE_DOSE",
            dose_mg=dose_mg,
            bioavailability_f=bioavailability_f,
            half_life_hr=round(t_half, 2),
            elimination_rate_ke=round(ke_hr, 4),
            absorption_rate_ka=round(ka_hr, 4),
            volume_distribution_l=round(vd_l, 2),
            clearance_l_hr=round(cl_l_hr, 2),
            cmax_mg_l=round(c_max, 4),
            tmax_hr=round(t_max, 2),
            auc_0_inf_mg_hr_l=round(auc_inf, 2),
            concentration_curve=curve,
        )

    @classmethod
    def simulate_iv_bolus(
        cls,
        dose_mg: float,
        ke_hr: float,
        vd_l: float,
        duration_hr: float = 24.0,
        num_points: int = 100,
    ) -> PKSimulationResult:
        if dose_mg <= 0 or ke_hr <= 0 or vd_l <= 0:
            raise ValueError("Dose, ke, and Vd must be positive for IV bolus simulation.")

        cl_l_hr = ke_hr * vd_l
        t_half = math.log(2.0) / ke_hr
        c0 = dose_mg / vd_l
        auc_inf = dose_mg / cl_l_hr

        dt = duration_hr / max(1, num_points - 1)
        curve = []
        for i in range(num_points):
            t = i * dt
            c = c0 * math.exp(-ke_hr * t)
            curve.append(PKSimulationPoint(time_hr=round(t, 2), plasma_conc_mg_l=round(c, 4)))

        return PKSimulationResult(
            dosing_route="IV_BOLUS",
            dose_mg=dose_mg,
            bioavailability_f=1.0,
            half_life_hr=round(t_half, 2),
            elimination_rate_ke=round(ke_hr, 4),
            absorption_rate_ka=None,
            volume_distribution_l=round(vd_l, 2),
            clearance_l_hr=round(cl_l_hr, 2),
            cmax_mg_l=round(c0, 4),
            tmax_hr=0.0,
            auc_0_inf_mg_hr_l=round(auc_inf, 2),
            concentration_curve=curve,
        )

    @classmethod
    def simulate_oral_multiple(
        cls,
        dose_mg: float,
        bioavailability_f: float,
        ka_hr: float,
        ke_hr: float,
        vd_l: float,
        dosing_interval_tau_hr: float = 12.0,
        num_doses: int = 7,
        num_points_per_interval: int = 25,
    ) -> PKSimulationResult:
        if dosing_interval_tau_hr <= 0 or num_doses <= 0:
            raise ValueError("Tau and num_doses must be positive.")
        if dose_mg <= 0 or bioavailability_f <= 0 or ka_hr <= 0 or ke_hr <= 0 or vd_l <= 0:
            raise ValueError("All PK parameters (dose, F, ka, ke, Vd) must be strictly positive.")

        # Handle near-equal ka and ke to avoid division by zero in Bateman function
        if abs(ka_hr - ke_hr) < 1e-6:
            ka_hr += 1e-4

        cl_l_hr = ke_hr * vd_l
        t_half = math.log(2.0) / ke_hr
        r_acc = 1.0 / (1.0 - math.exp(-ke_hr * dosing_interval_tau_hr))
        c_ss_avg = (dose_mg * bioavailability_f) / (cl_l_hr * dosing_interval_tau_hr)

        total_time = num_doses * dosing_interval_tau_hr
        dt = dosing_interval_tau_hr / num_points_per_interval
        total_points = int(total_time / dt) + 1

        curve = []
        max_seen = 0.0
        for p in range(total_points):
            t = p * dt
            c_tot = 0.0
            for d in range(num_doses):
                dose_time = d * dosing_interval_tau_hr
                if t >= dose_time:
                    t_since_dose = t - dose_time
                    c_dose = (dose_mg * bioavailability_f * ka_hr / (vd_l * (ka_hr - ke_hr))) * (
                        math.exp(-ke_hr * t_since_dose) - math.exp(-ka_hr * t_since_dose)
                    )
                    c_tot += max(0.0, c_dose)
            curve.append(PKSimulationPoint(time_hr=round(t, 2), plasma_conc_mg_l=round(c_tot, 4)))
            if c_tot > max_seen:
                max_seen = c_tot

        last_interval_pts = [pt.plasma_conc_mg_l for pt in curve if pt.time_hr >= (num_doses - 1) * dosing_interval_tau_hr]
        c_ss_max = max(last_interval_pts) if last_interval_pts else max_seen
        c_ss_min = min(last_interval_pts) if last_interval_pts else 0.0
        auc_tau = c_ss_avg * dosing_interval_tau_hr

        return PKSimulationResult(
            dosing_route="ORAL_MULTIPLE_DOSE",
            dose_mg=dose_mg,
            bioavailability_f=bioavailability_f,
            half_life_hr=round(t_half, 2),
            elimination_rate_ke=round(ke_hr, 4),
            absorption_rate_ka=round(ka_hr, 4),
            volume_distribution_l=round(vd_l, 2),
            clearance_l_hr=round(cl_l_hr, 2),
            cmax_mg_l=round(max_seen, 4),
            tmax_hr=round(math.log(ka_hr / ke_hr) / (ka_hr - ke_hr), 2),
            auc_0_inf_mg_hr_l=round((dose_mg * bioavailability_f) / cl_l_hr, 2),
            auc_tau_mg_hr_l=round(auc_tau, 2),
            c_ss_avg_mg_l=round(c_ss_avg, 4),
            c_ss_min_mg_l=round(c_ss_min, 4),
            c_ss_max_mg_l=round(c_ss_max, 4),
            accumulation_ratio=round(r_acc, 2),
            concentration_curve=curve,
        )


# ============================================================================
# Standard Reference Drug Database
# ============================================================================

REFERENCE_DRUGS = {
    "Aspirin": MoleculeProperties(
        name="Aspirin (Acetylsalicylic Acid)",
        mw=180.16, logp=1.19, hbd=1, hba=3, tpsa=63.6, rotatable_bonds=3,
        aromatic_rings=1, heavy_atoms=13, molar_refractivity=43.8, pka_acid=3.5, fsp3=0.11
    ),
    "Caffeine": MoleculeProperties(
        name="Caffeine",
        mw=194.19, logp=-0.07, hbd=0, hba=6, tpsa=58.4, rotatable_bonds=0,
        aromatic_rings=2, heavy_atoms=14, molar_refractivity=48.2, pka_base=0.6, fsp3=0.38
    ),
    "Ibuprofen": MoleculeProperties(
        name="Ibuprofen",
        mw=206.28, logp=3.50, hbd=1, hba=2, tpsa=37.3, rotatable_bonds=4,
        aromatic_rings=1, heavy_atoms=15, molar_refractivity=60.8, pka_acid=4.4, fsp3=0.54
    ),
    "Atorvastatin": MoleculeProperties(
        name="Atorvastatin",
        mw=558.64, logp=5.70, hbd=4, hba=7, tpsa=111.8, rotatable_bonds=12,
        aromatic_rings=4, heavy_atoms=41, molar_refractivity=159.2, pka_acid=4.5, fsp3=0.33
    ),
    "Imatinib": MoleculeProperties(
        name="Imatinib",
        mw=493.60, logp=3.50, hbd=2, hba=7, tpsa=86.3, rotatable_bonds=7,
        aromatic_rings=4, heavy_atoms=37, molar_refractivity=147.0, pka_base=8.1, fsp3=0.21
    ),
    "Morphine": MoleculeProperties(
        name="Morphine",
        mw=285.34, logp=0.89, hbd=2, hba=4, tpsa=49.8, rotatable_bonds=0,
        aromatic_rings=1, heavy_atoms=21, molar_refractivity=78.2, pka_base=8.2, fsp3=0.59
    ),
    "Vancomycin": MoleculeProperties(
        name="Vancomycin",
        mw=1449.25, logp=-3.10, hbd=19, hba=26, tpsa=470.0, rotatable_bonds=16,
        aromatic_rings=5, heavy_atoms=102, molar_refractivity=345.0, pka_base=7.75, fsp3=0.45
    ),
    "Metformin": MoleculeProperties(
        name="Metformin",
        mw=129.16, logp=-1.43, hbd=4, hba=3, tpsa=88.0, rotatable_bonds=2,
        aromatic_rings=0, heavy_atoms=9, molar_refractivity=35.0, pka_base=12.4, fsp3=0.50
    ),
    "Diazepam": MoleculeProperties(
        name="Diazepam",
        mw=284.74, logp=2.82, hbd=0, hba=3, tpsa=32.7, rotatable_bonds=1,
        aromatic_rings=2, heavy_atoms=20, molar_refractivity=80.1, pka_base=3.4, fsp3=0.12
    ),
}
