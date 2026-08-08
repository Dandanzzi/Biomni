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

# 에이전트 모드에서는 사용자의 자연어 질문을 그대로 A1에 넘긴다. 어떤 도구를 어떤 순서로 쓸지는
# 에이전트가 스스로 계획해야 하며(= Biomni의 핵심 동작), 아래 힌트는 --guided 플래그로만 덧붙는다.
WORKFLOW_HINT = (
    "\n\n[참고] 합성치사 전용 도구가 등록되어 있습니다: discover_synthetic_lethal_candidates, "
    "validate_sl_candidates_with_pubmed, analyze_ppi_network_for_sl, check_dependency_confounders, "
    "generate_sl_evidence_dossier. 근거가 약한 후보를 확정된 것처럼 제시하지 말고, 반대 근거와 QC 경고를 "
    "명시적으로 함께 보고하세요."
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


def run_agent(llm: str, data_path: str, query: str | None, guided: bool) -> None:
    """자연어 질문을 그대로 A1에게 넘긴다. 도구 선택과 실행 순서는 에이전트가 스스로 계획한다.

    query가 없으면 대화형으로 질문을 입력받는다 (빈 줄 입력 시 종료).
    """
    from dotenv import load_dotenv

    load_dotenv(".env")
    from biomni.agent import A1

    if not any(os.getenv(key) for key in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "AZURE_OPENAI_API_KEY")):
        raise SystemExit("LLM API 키가 없습니다. ANTHROPIC_API_KEY 등을 먼저 설정하세요.")

    agent = A1(path=data_path, llm=llm)
    registered = [tool["name"] for tool in agent.module2api.get("biomni.tool.synthetic_lethality", [])]
    print(f"에이전트에 등록된 합성치사 도구: {registered}")
    print("Biomni의 tool retriever가 질문에 맞는 도구를 스스로 선택합니다.\n")

    interactive = query is None
    while True:
        if interactive:
            try:
                question = input("질문> ").strip()
            except EOFError:
                break
            if not question:
                break
        else:
            question = query

        _, answer = agent.go(question + (WORKFLOW_HINT if guided else ""))
        print("\n" + "=" * 78)
        print("최종 답변")
        print("=" * 78)
        print(answer)

        if not interactive:
            break


def main() -> None:
    parser = argparse.ArgumentParser(description="KRAS-mutant pancreatic cancer synthetic lethality pipeline")
    parser.add_argument("--mode", choices=["direct", "agent"], default="direct")
    parser.add_argument("--cancer-type", default="Pancreatic Cancer")
    parser.add_argument("--mutation", default="KRAS")
    parser.add_argument("--top-n", type=int, default=5, help="candidates carried into validation")
    parser.add_argument("--output", default=None, help="write the full dossier to this file")
    parser.add_argument("--llm", default="claude-sonnet-4-5-20250929", help="agent mode only")
    parser.add_argument("--data-path", default="./data", help="agent mode only")
    parser.add_argument(
        "--query",
        default=None,
        help="에이전트에게 던질 자연어 질문. 생략하면 대화형으로 입력받는다 (agent 모드 전용)",
    )
    parser.add_argument(
        "--guided",
        action="store_true",
        help="등록된 합성치사 도구 목록을 힌트로 덧붙인다 (기본값은 에이전트가 스스로 계획)",
    )
    args = parser.parse_args()

    if args.mode == "direct":
        run_direct(args.cancer_type, args.mutation, args.top_n, args.output)
    else:
        run_agent(args.llm, args.data_path, args.query, args.guided)


if __name__ == "__main__":
    main()
