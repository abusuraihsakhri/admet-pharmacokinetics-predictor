#!/usr/bin/env python3
"""
Command-Line Interface for ADMET & Pharmacokinetics Predictor
=============================================================
Provides interactive and non-interactive workflows for evaluating drug-likeness,
QED score, CNS MPO score, ADMET risk profiles, and PK simulations.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import Optional, List

from admet_predictor import (
    MoleculeProperties,
    ADMETPredictor,
    PharmacokineticSimulator,
    REFERENCE_DRUGS,
    LipinskiRuleOf5,
    VeberRule,
    QEDCalculator,
    CNSMPOPredictor,
)


def format_report_text(report) -> str:
    """Format comprehensive ADMET report into a readable ASCII report."""
    m = report.molecule
    lines = [
        "=" * 78,
        f" ADMET & PHARMACOKINETICS CANDIDATE REPORT: {m.name}",
        "=" * 78,
        f"Physicochemical Properties:",
        f"  - Molecular Weight (MW):        {m.mw:.2f} Da",
        f"  - Octanol-Water Partition LogP: {m.logp:.2f}",
        f"  - H-Bond Donors (HBD):          {m.hbd}",
        f"  - H-Bond Acceptors (HBA):       {m.hba}",
        f"  - Polar Surface Area (TPSA):    {m.tpsa:.1f} Å²",
        f"  - Rotatable Bonds:              {m.rotatable_bonds}",
        f"  - Aromatic Rings:               {m.aromatic_rings}",
        f"  - Heavy Atoms:                  {m.heavy_atoms}",
        f"  - Molar Refractivity (MR):      {m.molar_refractivity:.1f}",
        f"  - Most Basic pKa:               {m.pka_base if m.pka_base is not None else 'N/A'}",
        f"  - Most Acidic pKa:              {m.pka_acid if m.pka_acid is not None else 'N/A'}",
        "-" * 78,
        "Drug-Likeness Rules & Filters:",
        f"  - Lipinski Rule of 5:  {'PASS' if report.lipinski.passes else 'FAIL'} ({report.lipinski.violations} violation(s)) - {report.lipinski.rationale}",
        f"  - Veber Filter:         {'PASS' if report.veber.passes else 'FAIL'} ({report.veber.violations} violation(s)) - {report.veber.rationale}",
        f"  - Egan Filter:          {'PASS' if report.egan.passes else 'FAIL'} - {report.egan.rationale}",
        f"  - Ghose Filter:         {'PASS' if report.ghose.passes else 'FAIL'} ({report.ghose.violations} violation(s))",
        f"  - Muegge (Bayer):       {'PASS' if report.muegge.passes else 'FAIL'} ({report.muegge.violations} violation(s))",
        f"  - Lead-Likeness:        {'PASS' if report.lead_likeness.passes else 'FAIL'}",
        "-" * 78,
        "Quantitative Metrics:",
        f"  - QED Score:            {report.qed.qed_score:.3f} [{report.qed.druglikeness_grade}]",
        f"  - CNS MPO Score:        {report.cns_mpo.score:.2f} / 6.00 [{report.cns_mpo.cns_permeability_likelihood}]",
        f"  - Predicted LogBB:      {report.cns_mpo.logbb_pred:.3f}",
        "-" * 78,
        "ADMET Profiling Estimates:",
        f"  - Human Intestinal Absorption (HIA):  {report.admet_prediction.hia_pct:.1f}%",
        f"  - Caco-2 Apparent Permeability:       {report.admet_prediction.caco2_perm_cm_s:.2f} x 10^-6 cm/s ({report.admet_prediction.caco2_class})",
        f"  - Plasma Protein Binding (PPB):       {report.admet_prediction.ppb_pct:.1f}% (Fraction unbound fu = {report.admet_prediction.fraction_unbound:.3f})",
        f"  - Volume of Distribution (Vd_ss):     {report.admet_prediction.vd_ss_l_kg:.2f} L/kg ({report.admet_prediction.vd_classification})",
        f"  - Total Clearance (CL):               {report.admet_prediction.clearance_ml_min_kg:.2f} mL/min/kg",
        f"  - Elimination Half-Life (t1/2):       {report.admet_prediction.elimination_half_life_hr:.2f} hours",
        f"  - CYP450 Inhibitions:                 {report.admet_prediction.cyp_inhibitions}",
        f"  - P-gp Substrate / Inhibitor:         {report.admet_prediction.pgp_substrate_risk} / {report.admet_prediction.pgp_inhibitor_risk}",
        f"  - Safety: hERG Cardiotoxicity:        {report.admet_prediction.herg_risk}",
        f"  - Safety: DILI Hepatotoxicity:        {report.admet_prediction.dili_risk}",
        f"  - Safety: Ames Mutagenicity:          {report.admet_prediction.ames_mutagenicity}",
        "=" * 78,
        f"OVERALL DEVELOPABILITY SCORE: {report.overall_druglikeness_score:.1f} / 100",
        f"Assessment: {report.overall_assessment}",
    ]
    if report.recommendations:
        lines.append("Actionable Recommendations:")
        for r in report.recommendations:
            lines.append(f"  * {r}")
    lines.append("=" * 78)
    return "\n".join(lines)


def run_interactive():
    """Interactive CLI wizard for molecular analysis and PK simulation."""
    print("=" * 70)
    print(" ADMET & Pharmacokinetics Predictor - Interactive Studio")
    print("=" * 70)
    print("1. Evaluate a Reference Drug (Aspirin, Caffeine, Atorvastatin, Imatinib, etc.)")
    print("2. Enter Custom Molecular Descriptors for Full ADMET Evaluation")
    print("3. Run 1-Compartment Oral Pharmacokinetic Simulation")
    print("4. Run IV Bolus Pharmacokinetic Simulation")
    print("5. Run Multiple-Dose Steady-State PK Simulation")
    print("6. Run Full Reference Benchmark Demo")
    print("q. Exit")
    print("-" * 70)

    choice = input("Select an option [1-6, q]: ").strip()
    if choice in ("q", "quit", "exit"):
        return 0

    if choice == "1":
        print("\nAvailable reference drugs:")
        keys = list(REFERENCE_DRUGS.keys())
        for idx, k in enumerate(keys, 1):
            print(f"  {idx}. {k}")
        d_choice = input(f"Choose drug [1-{len(keys)}]: ").strip()
        try:
            d_idx = int(d_choice) - 1
            if 0 <= d_idx < len(keys):
                drug_mol = REFERENCE_DRUGS[keys[d_idx]]
                rep = ADMETPredictor.evaluate_candidate(drug_mol)
                print(format_report_text(rep))
            else:
                print("Invalid selection.")
        except ValueError:
            print("Invalid input.")

    elif choice == "2":
        try:
            name = input("Molecule Name [Candidate-X]: ").strip() or "Candidate-X"
            mw = float(input("Molecular Weight (MW, Da) [350.0]: ").strip() or "350.0")
            logp = float(input("LogP [2.5]: ").strip() or "2.5")
            hbd = int(input("H-Bond Donors (HBD) [2]: ").strip() or "2")
            hba = int(input("H-Bond Acceptors (HBA) [4]: ").strip() or "4")
            tpsa = float(input("TPSA (Å²) [65.0]: ").strip() or "65.0")
            rot = int(input("Rotatable Bonds [4]: ").strip() or "4")
            arom = int(input("Aromatic Rings [2]: ").strip() or "2")
            heavy = int(input("Heavy Atoms [24]: ").strip() or "24")
            mr = float(input("Molar Refractivity [80.0]: ").strip() or "80.0")
            pka_b_str = input("Most Basic pKa (optional, Enter to skip): ").strip()
            pka_b = float(pka_b_str) if pka_b_str else None

            mol = MoleculeProperties(
                name=name, mw=mw, logp=logp, hbd=hbd, hba=hba,
                tpsa=tpsa, rotatable_bonds=rot, aromatic_rings=arom,
                heavy_atoms=heavy, molar_refractivity=mr, pka_base=pka_b
            )
            rep = ADMETPredictor.evaluate_candidate(mol)
            print(format_report_text(rep))
        except Exception as e:
            print(f"Error: {e}")

    elif choice == "3":
        try:
            dose = float(input("Dose (mg) [100.0]: ").strip() or "100.0")
            f = float(input("Bioavailability F (0-1) [0.8]: ").strip() or "0.8")
            ka = float(input("Absorption rate ka (1/hr) [1.2]: ").strip() or "1.2")
            ke = float(input("Elimination rate ke (1/hr) [0.15]: ").strip() or "0.15")
            vd = float(input("Volume of distribution Vd (L) [20.0]: ").strip() or "20.0")
            dur = float(input("Simulation Duration (hr) [24.0]: ").strip() or "24.0")

            sim = PharmacokineticSimulator.simulate_oral_single(dose, f, ka, ke, vd, dur)
            print("\n" + "=" * 60)
            print(f" ORAL SINGLE-DOSE PK SIMULATION (Dose={dose}mg, F={f})")
            print("=" * 60)
            print(f"  - Elimination Half-life (t1/2): {sim.half_life_hr:.2f} hr")
            print(f"  - Time to Peak (Tmax):          {sim.tmax_hr:.2f} hr")
            print(f"  - Peak Concentration (Cmax):    {sim.cmax_mg_l:.4f} mg/L")
            print(f"  - Total Clearance (CL):         {sim.clearance_l_hr:.2f} L/hr")
            print(f"  - Area Under Curve (AUC_0_inf): {sim.auc_0_inf_mg_hr_l:.2f} mg*hr/L")
            print("\nConcentration-Time Profile (Sampled):")
            for pt in sim.concentration_curve[:: max(1, len(sim.concentration_curve) // 10)]:
                bar = "#" * int(pt.plasma_conc_mg_l * 10 / max(1e-3, sim.cmax_mg_l))
                print(f"  t = {pt.time_hr:5.1f} hr: {pt.plasma_conc_mg_l:6.3f} mg/L | {bar}")
            print("=" * 60)
        except Exception as e:
            print(f"Error: {e}")

    elif choice in ("4", "5", "6"):
        run_demo()

    return 0


def run_demo(as_json: bool = False):
    """Run full benchmark across reference drug database."""
    reports = []
    for name, mol in REFERENCE_DRUGS.items():
        rep = ADMETPredictor.evaluate_candidate(mol)
        reports.append(rep)

    if as_json:
        print(json.dumps([r.to_dict() for r in reports], indent=2))
        return 0

    print("=" * 100)
    print(" ADMET & PHARMACOKINETICS REFERENCE BENCHMARK SUITE")
    print("=" * 100)
    header = f"{'Drug Name':<28} | {'MW':<6} | {'LogP':<5} | {'Ro5':<4} | {'Veber':<5} | {'QED':<5} | {'CNS MPO':<7} | {'HIA %':<5} | {'PPB %':<5} | {'Score':<5}"
    print(header)
    print("-" * 100)
    for r in reports:
        m = r.molecule
        ro5_str = "PASS" if r.lipinski.passes else "FAIL"
        veber_str = "PASS" if r.veber.passes else "FAIL"
        row = f"{m.name[:28]:<28} | {m.mw:6.1f} | {m.logp:5.2f} | {ro5_str:<4} | {veber_str:<5} | {r.qed.qed_score:5.3f} | {r.cns_mpo.score:7.2f} | {r.admet_prediction.hia_pct:5.1f} | {r.admet_prediction.ppb_pct:5.1f} | {r.overall_druglikeness_score:5.1f}"
        print(row)
    print("=" * 100)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="admet_predictor",
        description="ADMET & Pharmacokinetics Predictor Engine",
    )
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive studio")
    parser.add_argument("--demo", action="store_true", help="Run reference drug library benchmark")
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    sub = parser.add_subparsers(dest="command", help="Subcommands")

    # Evaluate single molecule
    eval_p = sub.add_parser("evaluate", help="Evaluate a single molecule")
    eval_p.add_argument("--name", default="Candidate", help="Molecule name")
    eval_p.add_argument("--mw", type=float, required=True, help="Molecular Weight (Da)")
    eval_p.add_argument("--logp", type=float, required=True, help="LogP (octanol-water)")
    eval_p.add_argument("--hbd", type=int, default=2, help="Hydrogen bond donors")
    eval_p.add_argument("--hba", type=int, default=4, help="Hydrogen bond acceptors")
    eval_p.add_argument("--tpsa", type=float, default=60.0, help="Topological polar surface area (Å²)")
    eval_p.add_argument("--rotatable", type=int, default=4, help="Rotatable bonds")
    eval_p.add_argument("--aromatic", type=int, default=2, help="Aromatic rings")
    eval_p.add_argument("--heavy", type=int, default=22, help="Heavy atoms count")
    eval_p.add_argument("--mr", type=float, default=75.0, help="Molar refractivity")
    eval_p.add_argument("--pka-base", type=float, default=None, help="Most basic pKa")
    eval_p.add_argument("--pka-acid", type=float, default=None, help="Most acidic pKa")
    eval_p.add_argument("--json", action="store_true", help="Output results in JSON format")

    # Reference drug lookup
    ref_p = sub.add_parser("ref", help="Evaluate a reference drug by name")
    ref_p.add_argument("name", choices=list(REFERENCE_DRUGS.keys()), help="Reference drug name")
    ref_p.add_argument("--json", action="store_true", help="Output results in JSON format")

    # Batch CSV processing
    batch_p = sub.add_parser("batch", help="Batch evaluate molecules from a CSV file")
    batch_p.add_argument("--input", "-in", required=True, help="Input CSV path")
    batch_p.add_argument("--output", "-out", required=True, help="Output CSV path")

    # PK Simulation
    pk_p = sub.add_parser("pk-sim", help="Pharmacokinetic concentration-time simulation")
    pk_p.add_argument("--route", choices=["oral", "iv", "multi"], default="oral", help="Dosing route")
    pk_p.add_argument("--dose", type=float, default=100.0, help="Dose in mg")
    pk_p.add_argument("--f", type=float, default=0.85, help="Oral bioavailability fraction (0-1)")
    pk_p.add_argument("--ka", type=float, default=1.2, help="Absorption rate constant ka (1/hr)")
    pk_p.add_argument("--ke", type=float, default=0.15, help="Elimination rate constant ke (1/hr)")
    pk_p.add_argument("--vd", type=float, default=25.0, help="Volume of distribution Vd (L)")
    pk_p.add_argument("--tau", type=float, default=12.0, help="Dosing interval tau (hr) for multi-dose")
    pk_p.add_argument("--doses", type=int, default=7, help="Number of doses for multi-dose")
    pk_p.add_argument("--duration", type=float, default=24.0, help="Simulation duration (hr)")
    pk_p.add_argument("--json", action="store_true", help="Output results in JSON format")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.interactive or (not args.command and not args.demo):
        if not args.demo and (argv is None or len(argv) == 0):
            return run_interactive()

    if args.demo:
        return run_demo(as_json=args.json)

    if args.command == "ref":
        mol = REFERENCE_DRUGS[args.name]
        rep = ADMETPredictor.evaluate_candidate(mol)
        if args.json:
            print(rep.to_json())
        else:
            print(format_report_text(rep))
        return 0

    if args.command == "evaluate":
        mol = MoleculeProperties(
            name=args.name,
            mw=args.mw,
            logp=args.logp,
            hbd=args.hbd,
            hba=args.hba,
            tpsa=args.tpsa,
            rotatable_bonds=args.rotatable,
            aromatic_rings=args.aromatic,
            heavy_atoms=args.heavy,
            molar_refractivity=args.mr,
            pka_base=args.pka_base,
            pka_acid=args.pka_acid,
        )
        rep = ADMETPredictor.evaluate_candidate(mol)
        if args.json:
            print(rep.to_json())
        else:
            print(format_report_text(rep))
        return 0

    if args.command == "batch":
        try:
            with open(args.input, "r", newline="", encoding="utf-8-sig") as f_in:
                reader = csv.DictReader(f_in)
                rows = list(reader)
            out_rows = []
            for r in rows:
                mol = MoleculeProperties(
                    name=r.get("name", "Unknown"),
                    mw=float(r.get("mw", 300.0)),
                    logp=float(r.get("logp", 2.0)),
                    hbd=int(r.get("hbd", 2)),
                    hba=int(r.get("hba", 4)),
                    tpsa=float(r.get("tpsa", 60.0)),
                    rotatable_bonds=int(r.get("rotatable_bonds", 4)),
                    aromatic_rings=int(r.get("aromatic_rings", 2)),
                    heavy_atoms=int(r.get("heavy_atoms", 22)),
                    molar_refractivity=float(r.get("molar_refractivity", 70.0)),
                )
                rep = ADMETPredictor.evaluate_candidate(mol)
                out_rows.append({
                    "name": mol.name,
                    "mw": mol.mw,
                    "logp": mol.logp,
                    "lipinski_passes": rep.lipinski.passes,
                    "lipinski_violations": rep.lipinski.violations,
                    "veber_passes": rep.veber.passes,
                    "qed_score": rep.qed.qed_score,
                    "cns_mpo_score": rep.cns_mpo.score,
                    "hia_pct": rep.admet_prediction.hia_pct,
                    "ppb_pct": rep.admet_prediction.ppb_pct,
                    "vd_ss_l_kg": rep.admet_prediction.vd_ss_l_kg,
                    "clearance_ml_min_kg": rep.admet_prediction.clearance_ml_min_kg,
                    "half_life_hr": rep.admet_prediction.elimination_half_life_hr,
                    "herg_risk": rep.admet_prediction.herg_risk,
                    "overall_score": rep.overall_druglikeness_score,
                })
            with open(args.output, "w", newline="", encoding="utf-8") as f_out:
                if out_rows:
                    writer = csv.DictWriter(f_out, fieldnames=list(out_rows[0].keys()))
                    writer.writeheader()
                    writer.writerows(out_rows)
            print(f"Successfully processed {len(out_rows)} molecules to {args.output}")
            return 0
        except Exception as e:
            print(f"Batch processing error: {e}", file=sys.stderr)
            return 1

    if args.command == "pk-sim":
        if args.route == "oral":
            res = PharmacokineticSimulator.simulate_oral_single(
                dose_mg=args.dose,
                bioavailability_f=args.f,
                ka_hr=args.ka,
                ke_hr=args.ke,
                vd_l=args.vd,
                duration_hr=args.duration,
            )
        elif args.route == "iv":
            res = PharmacokineticSimulator.simulate_iv_bolus(
                dose_mg=args.dose,
                ke_hr=args.ke,
                vd_l=args.vd,
                duration_hr=args.duration,
            )
        else: # multi
            res = PharmacokineticSimulator.simulate_oral_multiple(
                dose_mg=args.dose,
                bioavailability_f=args.f,
                ka_hr=args.ka,
                ke_hr=args.ke,
                vd_l=args.vd,
                dosing_interval_tau_hr=args.tau,
                num_doses=args.doses,
            )

        if args.json:
            print(json.dumps(asdict(res), indent=2))
        else:
            print("=" * 60)
            print(f" PK SIMULATION RESULT [{res.dosing_route}]")
            print("=" * 60)
            print(f"  - Dose:                       {res.dose_mg:.1f} mg (F = {res.bioavailability_f:.2f})")
            print(f"  - Half-Life (t1/2):           {res.half_life_hr:.2f} hr")
            print(f"  - Volume of Distribution:     {res.volume_distribution_l:.2f} L")
            print(f"  - Clearance (CL):             {res.clearance_l_hr:.2f} L/hr")
            print(f"  - Peak Concentration (Cmax):  {res.cmax_mg_l:.4f} mg/L")
            print(f"  - Time to Peak (Tmax):        {res.tmax_hr:.2f} hr")
            print(f"  - Total AUC (0-inf):          {res.auc_0_inf_mg_hr_l:.2f} mg*hr/L")
            if res.c_ss_avg_mg_l is not None:
                print(f"  - Steady-State Average (Css): {res.c_ss_avg_mg_l:.4f} mg/L")
                print(f"  - Accumulation Ratio (R):     {res.accumulation_ratio:.2f}")
            print("=" * 60)
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
