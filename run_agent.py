"""췌장암 KRAS 합성치사 후보 발굴 -> 실제 유전자명 확인 -> PubMed 검증 -> 실험 대상 선정.

이전 버전은 Gene_A / Gene_B / Gene_C 라는 가상의 이름과 손으로 적은 의존성 점수를 사용했기 때문에
"Gene_B"에 대응하는 실제 유전자가 존재하지 않았고, 그 결과로는 문헌 검증도 실험도 할 수 없었습니다.
이 스크립트는 동일한 질문을 실제 DepMap CRISPR 데이터로 다시 풀어서, 실험실에서 바로 주문할 수 있는
HUGO 유전자 심볼을 출력합니다.

    python run_agent.py                      # 췌장암 / KRAS (기본값)
    python run_agent.py --mutation TP53      # 다른 driver 변이
    python run_agent.py --skip-pubmed        # 통계 단계만 (오프라인)
"""

import argparse
import os

from biomni.tool.synthetic_lethality import (
    check_dependency_confounders,
    discover_synthetic_lethal_candidates,
    validate_sl_candidates_with_pubmed,
)

# 실험으로 넘길 후보를 고르는 기준. LLM 판단이 아니라 명시적 규칙으로 고정한다.
MAX_PAN_ESSENTIAL_PCT = 40.0  # 전체 세포주에서 이 비율 이상 필수면 치료 window 확보가 어렵다
MAX_Q_VALUE = 0.10  # 다중검정 보정 후에도 살아남아야 한다
MIN_MUTANT_DEPENDENT_PCT = 30.0  # 변이 세포주 중 실제로 의존하는 비율


def select_experiment_targets(table, top_k: int = 3):
    """통계 결과표에서 wet-lab 검증 우선순위를 규칙 기반으로 선별한다."""
    eligible = table[
        (table["pct_all_lines_dependent"] < MAX_PAN_ESSENTIAL_PCT)
        & (table["q_value"] < MAX_Q_VALUE)
        & (table["pct_mutant_dependent"] >= MIN_MUTANT_DEPENDENT_PCT)
    ].copy()
    # 선택도(변이 세포주 의존 비율 - 전체 세포주 의존 비율)가 클수록 genotype 특이적이다.
    eligible["selectivity_gap"] = eligible["pct_mutant_dependent"] - eligible["pct_all_lines_dependent"]
    eligible = eligible.sort_values(["selectivity_gap", "effect_difference"], ascending=[False, True])
    return eligible.head(top_k), len(eligible)


def parse_confounder_verdicts(report: str) -> dict:
    """교란요인 리포트에서 유전자별 판정(CONFOUNDED / DRIVER-CONSISTENT 등)을 추출한다."""
    verdicts = {}
    current = None
    for line in report.splitlines():
        if line.startswith("### "):
            current = line[4:].strip()
        elif current and line.strip().startswith("VERDICT:"):
            verdicts[current] = line.split("VERDICT:", 1)[1].strip()
            current = None
    return verdicts


def print_experiment_plan(gene: str, row, mutation: str) -> None:
    print("=" * 78)
    print(f"실험 계획: {mutation} 변이 의존적 {gene} 합성치사 검증")
    print("=" * 78)
    print(f"  가설: {mutation} 변이 췌장암 세포는 {gene} 결손에 선택적으로 취약하다.")
    print(
        f"  통계 근거: 변이 세포주 gene effect {row['mutant_mean_effect']:.3f} vs "
        f"야생형 {row['wildtype_mean_effect']:.3f} (차이 {row['effect_difference']:.3f}, "
        f"q={row['q_value']:.4f}, 변이 세포주의 {row['pct_mutant_dependent']:.0f}%가 의존)"
    )
    print("")
    print("  1) 세포주 패널 (DepMap에서 실제로 스크리닝된 계열)")
    print(f"     - {mutation} 변이군: MIA PaCa-2 (G12C), PANC-1 (G12D), AsPC-1 (G12D), CFPAC-1 (G12V)")
    print(f"     - {mutation} 야생형군: BxPC-3")
    print("     * 췌장암 KRAS 야생형 세포주는 DepMap 전체에서 4종뿐이므로, 야생형군만으로 대조하지 말고")
    print("       isogenic 계통(HPNE/HPDE + KRAS G12D 도입, 또는 KRAS degron 계열)을 함께 사용할 것.")
    print("")
    print("  2) 단일 perturbation")
    print(f"     - {gene}에 대해 서로 다른 표적 서열의 sgRNA 2종 + non-targeting 대조")
    print("     - 판독: 10일 CellTiter-Glo 생존율, caspase-3/7, colony formation")
    print(f"     - 반드시 western blot 또는 RT-qPCR로 {gene} 녹아웃 효율을 확인 (miss 시 위음성)")
    print("")
    print("  3) 특이성 확인 (on-target 검증)")
    print(f"     - sgRNA 저항성 {gene} cDNA 재발현 시 표현형이 회복되어야 한다 (rescue)")
    print(f"     - 회복되지 않으면 off-target 효과이므로 {gene}은 기각")
    print("")
    print("  4) 판정 기준")
    print("     - 변이군 대 야생형군 생존율 차이 >= 2배, p < 0.05 (two-way ANOVA의 genotype x KO 교호작용)")
    print("     - sgRNA 2종에서 모두 재현, rescue로 회복")
    print("     - 비형질전환 세포(HPNE 등)에서 동등한 의존성이 나오면 therapeutic window 없음 -> No-go")
    print("")


def main() -> None:
    parser = argparse.ArgumentParser(description="췌장암 합성치사 후보 발굴 및 실험 대상 선정")
    parser.add_argument("--cancer-type", default="Pancreatic Cancer")
    parser.add_argument("--mutation", default="KRAS")
    parser.add_argument("--top-n", type=int, default=15, help="통계 단계에서 보고할 후보 수")
    parser.add_argument("--skip-pubmed", action="store_true", help="문헌 검증 생략 (오프라인 실행)")
    parser.add_argument("--csv", default="sl_candidates.csv", help="후보 표를 저장할 경로")
    args = parser.parse_args()

    print("=" * 78)
    print(f"사용자 질문: {args.cancer_type}에서 {args.mutation} 변이 기반 합성치사 후보를 찾아줘")
    print("=" * 78)
    print("")

    # ---- 1단계: 실제 DepMap 데이터로 후보 발굴 -------------------------------------------
    report = discover_synthetic_lethal_candidates(
        cancer_type=args.cancer_type,
        target_mutation=args.mutation,
        top_n=args.top_n,
        output_csv_path=args.csv,
    )
    print(report)

    if report.startswith("FAILURE") or not os.path.exists(args.csv):
        print("\n후보 발굴에 실패했습니다. 위의 실패 사유를 확인하세요.")
        return

    import pandas as pd

    table = pd.read_csv(args.csv)
    if table.empty:
        print("\n통계 필터를 통과한 후보가 없습니다.")
        return

    # ---- 2단계: 실험 우선순위 선별 ---------------------------------------------------------
    targets, n_eligible = select_experiment_targets(table)
    print("")
    print("=" * 78)
    print("실험 우선순위 선별")
    print("=" * 78)
    print(
        f"  선별 기준: pan-essential < {MAX_PAN_ESSENTIAL_PCT:.0f}%, q < {MAX_Q_VALUE}, "
        f"변이 세포주 의존 비율 >= {MIN_MUTANT_DEPENDENT_PCT:.0f}%"
    )
    print(f"  전체 후보 {len(table)}개 중 {n_eligible}개 통과")
    if targets.empty:
        print("  기준을 통과한 후보가 없습니다. 통계 근거만으로 실험을 시작하지 마세요.")
        return

    print("")
    print(f"{'순위':<6}{'유전자':<12}{'변이군':>9}{'야생형':>9}{'차이':>9}{'q':>9}{'변이의존':>9}{'전체의존':>9}")
    print("-" * 78)
    for index, (_, row) in enumerate(targets.iterrows(), start=1):
        print(
            f"{index:<6}{row['gene']:<12}{row['mutant_mean_effect']:>9.3f}{row['wildtype_mean_effect']:>9.3f}"
            f"{row['effect_difference']:>9.3f}{row['q_value']:>9.4f}"
            f"{row['pct_mutant_dependent']:>8.0f}%{row['pct_all_lines_dependent']:>8.0f}%"
        )

    candidate_genes = targets["gene"].tolist()
    print("")
    print(f"=> 문헌 검증 및 실험 대상: {', '.join(candidate_genes)}")
    print("   (이전 버전의 'Gene_B'와 달리 모두 실존하는 HUGO 심볼이며 그대로 sgRNA 주문에 사용할 수 있습니다)")
    print("")

    # ---- 3단계: PubMed 검증 ----------------------------------------------------------------
    if not args.skip_pubmed:
        print(
            validate_sl_candidates_with_pubmed(
                disease=args.cancer_type,
                mutated_gene=args.mutation,
                candidate_genes=candidate_genes,
                max_papers_per_gene=6,
                include_abstracts=False,
            )
        )
        print("")

    # ---- 4단계: 교란요인 반증 검사 ---------------------------------------------------------
    # 실험비를 쓰기 전에, 이 의존성이 정말 driver 변이 때문인지 아니면 파라로그 손실이나
    # lineage 구성 때문인지 확인한다.
    confounder_report = check_dependency_confounders(
        target_mutation=args.mutation,
        candidate_genes=candidate_genes,
        cancer_type=args.cancer_type,
    )
    print(confounder_report)
    verdicts = parse_confounder_verdicts(confounder_report)

    # ---- 5단계: 최우선 후보의 실험 프로토콜 ------------------------------------------------
    clean = [g for g in candidate_genes if not verdicts.get(g, "").startswith("CONFOUNDED")]
    dropped = [g for g in candidate_genes if g not in clean]
    print("")
    if dropped:
        print("=" * 78)
        print("교란요인으로 실험 대상에서 제외된 후보")
        print("=" * 78)
        for gene in dropped:
            print(f"  {gene}: {verdicts[gene]}")
        print("")

    if not clean:
        print("모든 후보가 교란요인으로 설명됩니다. 실험을 시작하기 전에 해당 바이오마커로")
        print("층화한 뒤 discover_synthetic_lethal_candidates를 다시 실행하세요.")
        print(f"\n전체 후보 표: {args.csv}")
        return

    top = targets[targets["gene"] == clean[0]].iloc[0]
    print_experiment_plan(top["gene"], top, args.mutation)
    print(f"교란요인 판정: {verdicts.get(top['gene'], 'N/A')}")
    print(f"전체 후보 표: {args.csv}")


if __name__ == "__main__":
    main()
