"""Pancreatic cancer synthetic lethality scenario for the Biomni agent.

Two modes:

  direct  Run the tool chain deterministically, without an LLM. Requires only the local DepMap
          snapshot plus internet access to cBioPortal, NCBI Entrez and STRING.

              python run_sl_agent.py --mode direct

  agent   Hand the natural-language question to the Biomni A1 agent, which plans the tool calls
          itself. Requires an LLM API key (e.g. ANTHROPIC_API_KEY).

              python run_sl_agent.py --mode agent

The four tools live in biomni/tool/synthetic_lethality.py and are registered for the agent through
biomni/tool/tool_description/synthetic_lethality.py plus the field list in biomni/utils.py.
"""

import argparse
import os

USER_QUESTION = (
    "췌장암에서 KRAS 변이 기반 합성치사 후보를 찾고, 해당 후보들의 가능성을 논문으로 검증해 줘. "
    "발굴된 후보와 KRAS 사이의 단백질 상호작용 네트워크도 확인하고, "
    "통계적 근거 / 논문 근거 / 네트워크 근거를 통합한 종합 리포트를 만들어 줘."
)

AGENT_PROMPT = (
    "Find synthetic lethality candidates for KRAS-mutant pancreatic cancer and validate them.\n"
    "Follow this workflow:\n"
    "1. discover_synthetic_lethal_candidates(cancer_type='Pancreatic Cancer', target_mutation='KRAS') "
    "to obtain genotype-selective dependencies from DepMap CRISPR data.\n"
    "2. Take the genes on the CANDIDATE_GENES line (top 5) and run "
    "validate_sl_candidates_with_pubmed(disease='Pancreatic cancer', mutated_gene='KRAS', "
    "candidate_genes=<those genes>) for literature evidence.\n"
    "3. Run analyze_ppi_network_for_sl(target_gene='KRAS', candidate_genes=<those genes>) for "
    "protein-network evidence.\n"
    "4. Summarise as: discovered target -> statistical evidence -> literature evidence -> PPI evidence, "
    "and give a Go/Hold/No-go recommendation per candidate. State the contradicting evidence and the "
    "QC warnings explicitly - do not present weakly supported candidates as established.\n"
    "Alternatively, generate_sl_evidence_dossier(cancer_type='Pancreatic Cancer', target_mutation='KRAS') "
    "runs all three stages and integrates them in one call."
)


def run_direct(cancer_type: str, mutation: str, top_n: int, output_path: str | None) -> None:
    """Run the three-tool chain explicitly so each stage's raw output is visible."""
    from biomni.tool.synthetic_lethality import (
        analyze_ppi_network_for_sl,
        discover_synthetic_lethal_candidates,
        generate_sl_evidence_dossier,
        validate_sl_candidates_with_pubmed,
    )

    print(f"User question:\n  {USER_QUESTION}\n")
    print("=" * 78)
    print(f"STAGE 1/4 - DepMap statistical discovery ({mutation} in {cancer_type})")
    print("=" * 78)
    discovery = discover_synthetic_lethal_candidates(cancer_type=cancer_type, target_mutation=mutation, top_n=15)
    print(discovery)

    candidates = []
    for line in discovery.splitlines():
        if line.startswith("CANDIDATE_GENES:"):
            candidates = [g.strip() for g in line.split(":", 1)[1].split(",") if g.strip()][:top_n]
    if not candidates:
        print("\nNo candidate survived the statistical filters; stopping here.")
        return
    print(f"\n-> Carrying {len(candidates)} candidates into validation: {', '.join(candidates)}\n")

    print("=" * 78)
    print("STAGE 2/4 - PubMed literature validation (NCBI Entrez)")
    print("=" * 78)
    print(validate_sl_candidates_with_pubmed(cancer_type, mutation, candidates, max_papers_per_gene=6))

    print("\n" + "=" * 78)
    print("STAGE 3/4 - STRING protein-protein interaction analysis")
    print("=" * 78)
    print(analyze_ppi_network_for_sl(mutation, candidates))

    print("\n" + "=" * 78)
    print("STAGE 4/4 - Integrated evidence dossier")
    print("=" * 78)
    dossier = generate_sl_evidence_dossier(cancer_type=cancer_type, target_mutation=mutation, top_n=top_n)
    integrated = dossier.split("EVIDENCE DOSSIER")[-1]
    print("EVIDENCE DOSSIER" + integrated)

    if output_path:
        with open(output_path, "w") as handle:
            handle.write(dossier)
        print(f"\nFull dossier written to {output_path}")


def run_agent(llm: str, data_path: str) -> None:
    """Let the Biomni A1 agent plan and execute the tool chain from the natural-language question."""
    from biomni.agent import A1

    if not any(os.getenv(key) for key in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "AZURE_OPENAI_API_KEY")):
        raise SystemExit("No LLM API key found. Set ANTHROPIC_API_KEY (or another provider key) first.")

    agent = A1(path=data_path, llm=llm)
    registered = [tool["name"] for tool in agent.module2api.get("biomni.tool.synthetic_lethality", [])]
    print(f"Synthetic lethality tools available to the agent: {registered}\n")

    # go_stream already pretty-prints every step, so only the final answer is echoed here.
    final_output = None
    for step in agent.go_stream(AGENT_PROMPT):
        final_output = step.get("output")

    print("\n" + "=" * 78)
    print("FINAL AGENT ANSWER")
    print("=" * 78)
    print(final_output)


def main() -> None:
    parser = argparse.ArgumentParser(description="KRAS-mutant pancreatic cancer synthetic lethality pipeline")
    parser.add_argument("--mode", choices=["direct", "agent"], default="direct")
    parser.add_argument("--cancer-type", default="Pancreatic Cancer")
    parser.add_argument("--mutation", default="KRAS")
    parser.add_argument("--top-n", type=int, default=5, help="candidates carried into validation")
    parser.add_argument("--output", default=None, help="write the full dossier to this file")
    parser.add_argument("--llm", default="claude-sonnet-4-20250514", help="agent mode only")
    parser.add_argument("--data-path", default="./data", help="agent mode only")
    args = parser.parse_args()

    if args.mode == "direct":
        run_direct(args.cancer_type, args.mutation, args.top_n, args.output)
    else:
        run_agent(args.llm, args.data_path)


if __name__ == "__main__":
    main()
