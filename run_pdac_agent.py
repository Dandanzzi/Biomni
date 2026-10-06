"""PDAC translational pipeline: genomic data -> synthetic lethality -> drug -> organoid -> verdict.

Research question
-----------------
Can Biomni discover genetic vulnerabilities and synthetic lethal relationships in pancreatic cancer,
and can the predicted drug effects be functionally validated in patient-derived organoids (PDOs)?

The point of the chain is that the last stage can say NO. Stages 1-6 produce predictions; stages 7-8
read back what the wet lab measured and score the predictions against it, falsified ones first.

    1 driver landscape        profile_pdac_driver_landscape        which driver is even testable
    2 SL discovery            discover_synthetic_lethal_candidates  DepMap genotype contrast
    3 multi-channel ranking   discover_sl_multichannel              DAISY-style, rank not FDR gate
    4 organoid prioritisation rank_organoid_sl_candidates           which candidate to spend a PDO on
    5 drug mapping            map_sl_candidates_to_drugs            target -> compound, or CRISPR-only
    6 experiment design       design_organoid_sl_experiment         protocol with falsification criteria
    7 measured PDO response   analyze_organoid_drug_response        4PL fit, genotype comparison
      measured CRISPR         analyze_crispr_validation             knockout effect + selectivity
    8 scoring                 compare_prediction_with_experiment    precision@k, falsified list, verdict

Two modes:

  direct  Run the chain deterministically, without an LLM:

              python run_pdac_agent.py --mode direct --stages 1-6
              python run_pdac_agent.py --mode direct --stages 7-8 \
                  --screen-csv pdo_screen.csv --genotype-csv pdo_genotype.csv \
                  --crispr-csv crispr_validation.csv

  agent   Hand the research question to the Biomni A1 agent and let it plan the tool calls:

              python run_pdac_agent.py --mode agent
"""

import argparse
import os

RESEARCH_QUESTION = (
    "췌장암(PDAC)의 유전적 취약성과 합성치사 관계를 발굴하고, 예측된 약물의 치료 효과를 "
    "환자유래 오가노이드(PDO)에서 기능적으로 검증할 수 있는가? "
    "DepMap 코호트에서 어떤 드라이버가 실제로 검정 가능한지부터 확인하고, 합성치사 후보를 도출한 뒤, "
    "각 후보가 약물로 공략 가능한지 아니면 CRISPR 전용인지 구분하고, 오가노이드 검증 실험을 설계해 줘. "
    "근거가 약한 후보를 확정된 것처럼 제시하지 말고 반증 근거와 QC 경고를 함께 보고해."
)

WORKFLOW_HINT = (
    "\n\n[참고] PDAC 변환연구 도구가 등록되어 있습니다: profile_pdac_driver_landscape, "
    "map_sl_candidates_to_drugs, analyze_organoid_drug_response, analyze_crispr_validation, "
    "compare_prediction_with_experiment. 합성치사 발굴/오가노이드 도구와 함께 사용하세요."
)

DEFAULT_OUTPUT_DIR = "./pdac_run"


def _section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def _parse_candidates(report: str, marker: str = "CANDIDATE_GENES:") -> list[str]:
    """Pull a machine-readable gene list out of a tool report."""
    for line in report.splitlines():
        if line.strip().startswith(marker):
            payload = line.split(":", 1)[1]
            return [gene.strip() for gene in payload.split(",") if gene.strip()]
    return []


def run_prediction_stages(args) -> None:
    """Stages 1-6: everything Biomni can do before the lab is involved."""
    from biomni.tool.organoid_sl import design_organoid_sl_experiment, rank_organoid_sl_candidates
    from biomni.tool.pdac_translation import map_sl_candidates_to_drugs, profile_pdac_driver_landscape
    from biomni.tool.synthetic_lethality import discover_synthetic_lethal_candidates

    os.makedirs(args.output_dir, exist_ok=True)
    prediction_csv = os.path.join(args.output_dir, "prediction_depmap.csv")

    _section(f"STAGE 1/6 - PDAC driver landscape (which driver can this cohort test?)")
    print(
        profile_pdac_driver_landscape(
            cancer_type=args.cancer_type,
            output_csv_path=os.path.join(args.output_dir, "driver_landscape.csv"),
        )
    )

    _section(f"STAGE 2/6 - Synthetic lethal discovery ({args.driver} in {args.cancer_type})")
    discovery = discover_synthetic_lethal_candidates(
        cancer_type=args.cancer_type,
        target_mutation=args.driver,
        top_n=args.top_n,
        output_csv_path=prediction_csv,
    )
    print(discovery)
    candidates = _parse_candidates(discovery)
    if not candidates:
        print(
            "\nNo candidate cleared the FDR gate. This is expected for KRAS in PDAC - the wild-type arm is "
            "tiny (see stage 1) - so continue with the rank-based channel instead of stopping here."
        )

    if args.multichannel:
        _section("STAGE 3/6 - Multi-channel ranking (DAISY-style, no FDR gate)")
        from biomni.tool.sl_multichannel import discover_sl_multichannel

        ranked = discover_sl_multichannel(
            driver_genes=[args.driver], cancer_type=args.cancer_type, top_k=args.top_n
        )
        print(ranked)
        candidates = candidates or _parse_candidates(ranked, "TOP_CANDIDATES:")
    else:
        print("\n(stage 3 skipped: pass --multichannel to run the rank-based channel)")

    if not candidates:
        candidates = [gene.strip().upper() for gene in args.fallback_candidates.split(",") if gene.strip()]
        print(f"\nUsing fallback candidate list for the downstream stages: {', '.join(candidates)}")

    _section("STAGE 4/6 - Organoid prioritisation (which candidate is worth a PDO?)")
    ranking = rank_organoid_sl_candidates(
        candidate_genes=candidates, target_mutation=args.driver, cancer_type=args.cancer_type, top_n=args.top_n
    )
    print(ranking)

    _section("STAGE 5/6 - Drug candidate identification")
    drugs = map_sl_candidates_to_drugs(
        candidate_genes=candidates,
        driver_gene=args.driver,
        output_csv_path=os.path.join(args.output_dir, "drug_candidates.csv"),
    )
    print(drugs)

    top_pick = candidates[0]
    for line in ranking.splitlines():
        if line.strip().startswith("ORGANOID_TOP_PICK:"):
            top_pick = line.split(":", 1)[1].strip().split()[0]
            break

    _section(f"STAGE 6/6 - Organoid experiment design ({top_pick})")
    print(
        design_organoid_sl_experiment(
            candidate_gene=top_pick, target_mutation=args.driver, cancer_type=args.cancer_type
        )
    )

    print("\n" + "-" * 78)
    print(f"Predictions written to {prediction_csv}")
    print(
        "Run the organoid experiment, then score the predictions with:\n"
        f"  python run_pdac_agent.py --mode direct --stages 7-8 --screen-csv <pdo_screen.csv> "
        f"--crispr-csv <crispr_validation.csv> --genotype-csv <pdo_genotype.csv> "
        f"--prediction-csv {prediction_csv}"
    )


def run_experiment_stages(args) -> None:
    """Stages 7-8: read back what the lab measured and score the predictions against it."""
    from biomni.tool.pdac_translation import (
        analyze_crispr_validation,
        analyze_organoid_drug_response,
        compare_prediction_with_experiment,
    )

    os.makedirs(args.output_dir, exist_ok=True)
    crispr_summary = os.path.join(args.output_dir, "crispr_summary.csv")

    if args.screen_csv:
        _section("STAGE 7a - Measured organoid drug response")
        print(
            analyze_organoid_drug_response(
                screen_csv_path=args.screen_csv,
                genotype_csv_path=args.genotype_csv,
                output_csv_path=os.path.join(args.output_dir, "pdo_curves.csv"),
                plot_output_path=os.path.join(args.output_dir, "pdo_doseresponse.png"),
            )
        )
    else:
        print("(stage 7a skipped: no --screen-csv)")

    if args.crispr_csv:
        _section("STAGE 7b - Measured organoid CRISPR validation")
        print(
            analyze_crispr_validation(
                validation_csv_path=args.crispr_csv,
                genotype_csv_path=args.genotype_csv,
                output_csv_path=crispr_summary,
            )
        )
    else:
        print("(stage 7b skipped: no --crispr-csv)")

    if args.crispr_csv and args.prediction_csv:
        _section("STAGE 8 - Biomni prediction vs experimental result")
        print(
            compare_prediction_with_experiment(
                prediction_csv_path=args.prediction_csv,
                experiment_csv_path=crispr_summary,
                top_k=args.top_k,
                output_csv_path=os.path.join(args.output_dir, "concordance.csv"),
            )
        )
    else:
        print(
            "\n(stage 8 skipped: it needs both --prediction-csv from stage 2 and --crispr-csv measurements)"
        )


def run_agent(args) -> None:
    """Hand the research question to A1 and let the agent plan the tool calls itself."""
    from dotenv import load_dotenv

    load_dotenv(".env")
    from biomni.agent import A1

    if not any(os.getenv(key) for key in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "AZURE_OPENAI_API_KEY")):
        raise SystemExit("LLM API 키가 없습니다. ANTHROPIC_API_KEY 등을 먼저 설정하세요.")

    agent = A1(path=args.data_path, llm=args.llm)
    for module in ("biomni.tool.pdac_translation", "biomni.tool.synthetic_lethality", "biomni.tool.organoid_sl"):
        registered = [tool["name"] for tool in agent.module2api.get(module, [])]
        print(f"{module}: {registered}")
    print()

    question = args.query or RESEARCH_QUESTION
    _, answer = agent.go(question + (WORKFLOW_HINT if args.guided else ""))
    _section("최종 답변")
    print(answer)


def main() -> None:
    parser = argparse.ArgumentParser(description="PDAC synthetic lethality -> organoid validation pipeline")
    parser.add_argument("--mode", choices=["direct", "agent"], default="direct")
    parser.add_argument("--stages", choices=["1-6", "7-8", "all"], default="1-6")
    parser.add_argument("--cancer-type", default="Pancreatic Cancer")
    parser.add_argument("--driver", default="KRAS")
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("--top-k", type=int, default=5, help="stage 8 precision@k")
    parser.add_argument("--multichannel", action="store_true", help="run the rank-based channel (stage 3)")
    parser.add_argument(
        "--fallback-candidates",
        default="TEAD1,SHOC2,PTPN11,WRN,CDK4",
        help="candidates used downstream when the FDR gate returns nothing",
    )
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--screen-csv", default=None, help="measured PDO drug screen (stage 7a)")
    parser.add_argument("--crispr-csv", default=None, help="measured organoid CRISPR validation (stage 7b)")
    parser.add_argument("--genotype-csv", default=None, help="organoid -> genotype mapping")
    parser.add_argument("--prediction-csv", default=None, help="prediction table from stage 2 (stage 8)")
    parser.add_argument("--llm", default="claude-sonnet-4-5-20250929", help="agent mode only")
    parser.add_argument("--data-path", default="./data", help="agent mode only")
    parser.add_argument("--query", default=None, help="override the research question (agent mode)")
    parser.add_argument("--guided", action="store_true", help="append the registered tool list as a hint")
    args = parser.parse_args()

    if args.mode == "agent":
        run_agent(args)
        return
    if args.stages in ("1-6", "all"):
        run_prediction_stages(args)
    if args.stages in ("7-8", "all"):
        run_experiment_stages(args)


if __name__ == "__main__":
    main()
