#!/usr/bin/env python3
"""
Unit Test Suite for ADMET & Pharmacokinetics Predictor
======================================================
Comprehensive test suite verifying physicochemical evaluations,
drug-likeness filters, QED, CNS MPO, ADMET properties, and PK models.
"""

import math
import os
import sys
import csv
import tempfile
import unittest
import json
from pathlib import Path

# Ensure root directory is on path
ROOT_DIR = Path(__file__).resolve().parent
if ROOT_DIR.name == "tests":
    ROOT_DIR = ROOT_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from admet_predictor import (
    MoleculeProperties,
    FilterResult,
    QEDResult,
    CNSMPOResult,
    ADMETPropertyEstimate,
    PKSimulationResult,
    LipinskiRuleOf5,
    VeberRule,
    EganRule,
    GhoseFilter,
    MueggeFilter,
    LeadLikenessFilter,
    QEDCalculator,
    CNSMPOPredictor,
    ADMETPredictor,
    PharmacokineticSimulator,
    REFERENCE_DRUGS,
)
import cli


class TestMoleculeProperties(unittest.TestCase):
    """Test validation and default computations of MoleculeProperties."""

    def test_valid_molecule_initialization(self):
        mol = MoleculeProperties(
            name="TestMol", mw=250.0, logp=1.8, hbd=1, hba=3,
            tpsa=45.0, rotatable_bonds=2, aromatic_rings=1, heavy_atoms=18
        )
        self.assertEqual(mol.name, "TestMol")
        self.assertEqual(mol.mw, 250.0)
        self.assertEqual(mol.logp, 1.8)
        self.assertEqual(mol.logd74, 1.8)

    def test_invalid_negative_mw_raises_value_error(self):
        with self.assertRaises(ValueError):
            MoleculeProperties(mw=-10.0)

    def test_invalid_negative_descriptors_raise_value_error(self):
        with self.assertRaises(ValueError):
            MoleculeProperties(hbd=-1)
        with self.assertRaises(ValueError):
            MoleculeProperties(tpsa=-5.0)
        with self.assertRaises(ValueError):
            MoleculeProperties(rotatable_bonds=-2)

    def test_basic_drug_logd74_ionization(self):
        mol = MoleculeProperties(name="BasicDrug", mw=200.0, logp=3.0, pka_base=9.4)
        self.assertLess(mol.logd74, mol.logp)
        self.assertAlmostEqual(mol.logd74, 3.0 - math.log10(1.0 + 10 ** (9.4 - 7.4)), places=1)

    def test_acidic_drug_logd74_ionization(self):
        mol = MoleculeProperties(name="AcidDrug", mw=200.0, logp=3.0, pka_acid=4.4)
        self.assertLess(mol.logd74, mol.logp)


class TestDrugLikenessFilters(unittest.TestCase):
    """Test standard rule-based drug-likeness filters."""

    def test_lipinski_pass_zero_violations(self):
        mol = REFERENCE_DRUGS["Aspirin"]
        res = LipinskiRuleOf5.evaluate(mol)
        self.assertTrue(res.passes)
        self.assertEqual(res.violations, 0)

    def test_lipinski_pass_one_violation(self):
        mol = MoleculeProperties(mw=520.0, logp=3.0, hbd=2, hba=6)
        res = LipinskiRuleOf5.evaluate(mol)
        self.assertTrue(res.passes)
        self.assertEqual(res.violations, 1)

    def test_lipinski_fail_multiple_violations(self):
        vanc = REFERENCE_DRUGS["Vancomycin"]
        res = LipinskiRuleOf5.evaluate(vanc)
        self.assertFalse(res.passes)
        self.assertGreaterEqual(res.violations, 3)

    def test_lipinski_boundary_exact_thresholds(self):
        boundary_mol = MoleculeProperties(mw=500.0, logp=5.0, hbd=5, hba=10)
        res = LipinskiRuleOf5.evaluate(boundary_mol)
        self.assertTrue(res.passes)
        self.assertEqual(res.violations, 0)

    def test_veber_rules(self):
        caffeine = REFERENCE_DRUGS["Caffeine"]
        res = VeberRule.evaluate(caffeine)
        self.assertTrue(res.passes)
        self.assertEqual(res.violations, 0)

        bulky = MoleculeProperties(rotatable_bonds=14, tpsa=160.0)
        res_bulky = VeberRule.evaluate(bulky)
        self.assertFalse(res_bulky.passes)
        self.assertEqual(res_bulky.violations, 2)

    def test_egan_filter(self):
        mol_pass = MoleculeProperties(logp=2.5, tpsa=70.0)
        self.assertTrue(EganRule.evaluate(mol_pass).passes)

        mol_fail_tpsa = MoleculeProperties(logp=2.5, tpsa=145.0)
        self.assertFalse(EganRule.evaluate(mol_fail_tpsa).passes)

    def test_ghose_filter(self):
        # Imatinib has MW 493.6, MR 147 (fails slightly), let's test a compliant molecule
        compliant = MoleculeProperties(mw=320.0, logp=2.4, molar_refractivity=80.0, heavy_atoms=24)
        res = GhoseFilter.evaluate(compliant)
        self.assertTrue(res.passes)
        self.assertEqual(res.violations, 0)

    def test_muegge_filter(self):
        aspirin = REFERENCE_DRUGS["Aspirin"]
        res = MueggeFilter.evaluate(aspirin)
        self.assertEqual(res.violations, 1)

    def test_lead_likeness_filter(self):
        lead_mol = MoleculeProperties(mw=220.0, logp=2.1, rotatable_bonds=4)
        self.assertTrue(LeadLikenessFilter.evaluate(lead_mol).passes)

        heavy_mol = MoleculeProperties(mw=480.0, logp=4.5, rotatable_bonds=9)
        self.assertFalse(LeadLikenessFilter.evaluate(heavy_mol).passes)


class TestQEDCalculator(unittest.TestCase):
    """Test Quantitative Estimate of Drug-likeness (QED) scoring."""

    def test_qed_ideal_range(self):
        asp_qed = QEDCalculator.calculate(REFERENCE_DRUGS["Aspirin"])
        self.assertGreater(asp_qed.qed_score, 0.45)
        self.assertLessEqual(asp_qed.qed_score, 1.0)
        self.assertIn("mw", asp_qed.individual_desirabilities)

    def test_qed_vancomycin_low(self):
        vanc_qed = QEDCalculator.calculate(REFERENCE_DRUGS["Vancomycin"])
        self.assertLess(vanc_qed.qed_score, 0.35)
        self.assertEqual(vanc_qed.druglikeness_grade, "Low Drug-Likeness")


class TestCNSMPOPredictor(unittest.TestCase):
    """Test Pfizer CNS Multiparameter Optimization scoring & BBB."""

    def test_cns_mpo_diazepam_high_permeability(self):
        diaz = REFERENCE_DRUGS["Diazepam"]
        res = CNSMPOPredictor.calculate(diaz)
        self.assertGreaterEqual(res.score, 4.0)
        self.assertEqual(res.cns_permeability_likelihood, "High (CNS+)")
        self.assertGreater(res.logbb_pred, -0.5)

    def test_cns_mpo_vancomycin_non_permeant(self):
        vanc = REFERENCE_DRUGS["Vancomycin"]
        res = CNSMPOPredictor.calculate(vanc)
        self.assertLessEqual(res.score, 3.0)
        self.assertIn("Low", res.cns_permeability_likelihood)
        self.assertGreater(len(res.risk_factors), 0)

    def test_cns_mpo_tpsa_ramp(self):
        tpsa_low = CNSMPOPredictor._tpsa_score(15.0)
        tpsa_opt = CNSMPOPredictor._tpsa_score(60.0)
        tpsa_high = CNSMPOPredictor._tpsa_score(140.0)
        self.assertEqual(tpsa_low, 0.0)
        self.assertEqual(tpsa_opt, 1.0)
        self.assertEqual(tpsa_high, 0.0)


class TestADMETPredictor(unittest.TestCase):
    """Test ADMET profiling parameters."""

    def test_admet_prediction_aspirin(self):
        asp = REFERENCE_DRUGS["Aspirin"]
        pred = ADMETPredictor.predict_properties(asp)
        self.assertGreater(pred.hia_pct, 30.0)
        self.assertGreater(pred.caco2_perm_cm_s, 5.0)
        self.assertLess(pred.vd_ss_l_kg, 1.5)
        self.assertEqual(pred.herg_risk, "Low")

    def test_comprehensive_report_generation(self):
        mol = REFERENCE_DRUGS["Imatinib"]
        rep = ADMETPredictor.evaluate_candidate(mol)
        self.assertIsInstance(rep.to_dict(), dict)
        self.assertIsInstance(rep.to_json(), str)
        self.assertGreaterEqual(rep.overall_druglikeness_score, 0.0)
        self.assertLessEqual(rep.overall_druglikeness_score, 100.0)


class TestPharmacokineticSimulator(unittest.TestCase):
    """Test analytical PK simulation formulas."""

    def test_oral_single_dose_analytical_solution(self):
        dose = 100.0
        f = 0.9
        ka = 1.5
        ke = 0.1
        vd = 20.0
        sim = PharmacokineticSimulator.simulate_oral_single(dose, f, ka, ke, vd, duration_hr=24.0)

        expected_t_half = math.log(2.0) / ke
        self.assertAlmostEqual(sim.half_life_hr, round(expected_t_half, 2), places=2)

        expected_tmax = math.log(ka / ke) / (ka - ke)
        self.assertAlmostEqual(sim.tmax_hr, round(expected_tmax, 2), places=2)

        expected_auc = (dose * f) / (ke * vd)
        self.assertAlmostEqual(sim.auc_0_inf_mg_hr_l, round(expected_auc, 2), places=2)
        self.assertEqual(sim.concentration_curve[0].plasma_conc_mg_l, 0.0)

    def test_iv_bolus_simulation(self):
        dose = 200.0
        ke = 0.2
        vd = 40.0
        sim = PharmacokineticSimulator.simulate_iv_bolus(dose, ke, vd, duration_hr=12.0)

        self.assertAlmostEqual(sim.cmax_mg_l, 5.0, places=2)
        self.assertEqual(sim.tmax_hr, 0.0)
        self.assertAlmostEqual(sim.concentration_curve[0].plasma_conc_mg_l, 5.0, places=2)

    def test_oral_multiple_dosing_steady_state(self):
        dose = 250.0
        f = 0.8
        ka = 1.0
        ke = 0.1
        vd = 30.0
        tau = 12.0
        sim = PharmacokineticSimulator.simulate_oral_multiple(dose, f, ka, ke, vd, dosing_interval_tau_hr=tau, num_doses=6)

        expected_r = 1.0 / (1.0 - math.exp(-ke * tau))
        self.assertAlmostEqual(sim.accumulation_ratio, round(expected_r, 2), places=2)
        self.assertIsNotNone(sim.c_ss_avg_mg_l)
        self.assertGreater(sim.c_ss_max_mg_l, sim.c_ss_min_mg_l)

    def test_invalid_pk_parameters_raise_error(self):
        with self.assertRaises(ValueError):
            PharmacokineticSimulator.simulate_oral_single(dose_mg=-10, bioavailability_f=0.8, ka_hr=1, ke_hr=0.1, vd_l=10)


class TestCLIAndBatch(unittest.TestCase):
    """Test CLI workflows, JSON formatting, and batch CSV processing."""

    def test_cli_demo_execution(self):
        res = cli.main(["--demo"])
        self.assertEqual(res, 0)

    def test_cli_ref_lookup(self):
        res = cli.main(["ref", "Aspirin"])
        self.assertEqual(res, 0)

    def test_cli_evaluate_subcommand(self):
        res = cli.main(["evaluate", "--mw", "320", "--logp", "2.1", "--hbd", "1", "--hba", "4", "--name", "Candidate-101"])
        self.assertEqual(res, 0)

    def test_cli_pk_simulation(self):
        res = cli.main(["pk-sim", "--route", "oral", "--dose", "150", "--f", "0.85"])
        self.assertEqual(res, 0)

    def test_batch_csv_processing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_csv = os.path.join(tmpdir, "input.csv")
            out_csv = os.path.join(tmpdir, "output.csv")
            with open(in_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["name", "mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "aromatic_rings", "heavy_atoms", "molar_refractivity"])
                writer.writerow(["DrugA", 250.0, 1.8, 1, 3, 50.0, 3, 1, 18, 65.0])
                writer.writerow(["DrugB", 580.0, 5.5, 6, 12, 160.0, 12, 4, 45, 160.0])

            ret = cli.main(["batch", "--input", in_csv, "--output", out_csv])
            self.assertEqual(ret, 0)
            self.assertTrue(os.path.exists(out_csv))
            with open(out_csv, "r") as f_out:
                lines = f_out.readlines()
                self.assertEqual(len(lines), 3)


if __name__ == "__main__":
    unittest.main()
