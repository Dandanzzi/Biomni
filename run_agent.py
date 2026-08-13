"""췌장암 KRAS 합성치사: 세포주 관점과 오가노이드 관점을 나란히 산출한다.

LLM을 쓰지 않는 결정론적 파이프라인이라 같은 입력이면 항상 같은 결과가 나옵니다.
(에이전트가 스스로 도구를 계획하는 버전은 run_sl_agent.py --mode agent)

두 관점을 모두 출력합니다.

  [세포주 관점]  DepMap CRISPR 47개 췌장암 세포주에서 통계적으로 발굴한 후보.
                 근거는 가장 강하지만, 2D 배양·아형 평균화·정상 대조군 부재라는 한계를 안습니다.
  [오가노이드 관점] 세포주 패널이 원리적으로 답할 수 없는 후보 — 아형 평균화로 상쇄된 것들 —
                 그리고 각 후보가 오가노이드 배지에서 가려질지(위음성) 여부.

    python run_agent.py                      # 췌장암 / KRAS (기본값)
    python run_agent.py --mutation TP53      # 다른 driver 변이
    python run_agent.py --skip-pubmed        # 문헌 검증 생략
    python run_agent.py --skip-organoid      # 세포주 관점만
"""

import argparse
import os

from biomni.tool.organoid_sl import (
    assess_organoid_transferability,
    design_organoid_sl_experiment,
    discover_subtype_masked_sl_candidates,
)
from biomni.tool.synthetic_lethality import (
    check_dependency_confounders,
    discover_synthetic_lethal_candidates,
    validate_sl_candidates_with_pubmed,
)

# 실험으로 넘길 후보를 고르는 기준. LLM 판단이 아니라 명시적 규칙으로 고정한다.
MAX_PAN_ESSENTIAL_PCT = 40.0  # 전체 세포주에서 이 비율 이상 필수면 치료 window 확보가 어렵다
MAX_Q_VALUE = 0.10  # 다중검정 보정 후에도 살아남아야 한다
MIN_MUTANT_DEPENDENT_PCT = 30.0  # 변이 세포주 중 실제로 의존하는 비율


def print_interpretation_guide(table, targets) -> None:
    """후보표를 어떻게 읽어야 하는지, 왜 이 후보가 선택되었는지를 화면에 설명한다."""
    print("=" * 78)
    print("후보표 해석")
    print("=" * 78)
    print("  effect_difference  변이군 - 야생형군 gene effect. 음수일수록 합성치사 방향.")
    print("  q_value            17,787개 유전자 다중검정 보정값. p값이 아니라 이 값을 보세요.")
    print("  pct_wildtype_dep   야생형 세포주 중 의존 비율. 0%여야 genotype 선택적입니다.")
    print("  pct_all_dep        전체 1,183개 세포주 중 의존 비율. 높으면 그냥 필수유전자이고")
    print("                     정상세포도 죽으므로 치료 window가 없습니다.")
    print("  선택도             변이의존% - 전체의존%. 이 값이 실질적인 우선순위 지표입니다.")
    print("")
    print(f"{'유전자':<10}{'diff':>8}{'q':>9}{'%WT의존':>9}{'%전체의존':>10}{'선택도':>9}  판정")
    print("-" * 78)
    chosen = set(targets["gene"])
    for _, row in table.sort_values("effect_difference").iterrows():
        gap = row["pct_mutant_dependent"] - row["pct_all_lines_dependent"]
        if row["gene"] in chosen:
            reason = "선별 통과"
        elif row["q_value"] >= MAX_Q_VALUE:
            reason = f"FDR 탈락 (q={row['q_value']:.2f})"
        elif row["pct_all_lines_dependent"] >= MAX_PAN_ESSENTIAL_PCT:
            reason = f"전체 {row['pct_all_lines_dependent']:.0f}% 의존 - 치료 window 없음"
        elif row["pct_mutant_dependent"] < MIN_MUTANT_DEPENDENT_PCT:
            reason = f"변이주 중 {row['pct_mutant_dependent']:.0f}%만 의존 - 효과 미미"
        else:
            reason = "-"
        if row["pct_wildtype_dependent"] > 0:
            reason += f" / 야생형도 {row['pct_wildtype_dependent']:.0f}% 의존"
        print(
            f"{row['gene']:<10}{row['effect_difference']:>8.3f}{row['q_value']:>9.3f}"
            f"{row['pct_wildtype_dependent']:>8.0f}%{row['pct_all_lines_dependent']:>9.0f}%{gap:>9.1f}  {reason}"
        )
    print("")


def select_experiment_targets(table, top_k: int = 5):
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


def parse_labelled_block(report: str, label: str) -> dict:
    """'### GENE' 블록에서 'label: ...' 줄을 유전자별로 뽑아낸다."""
    found = {}
    current = None
    for line in report.splitlines():
        if line.startswith("### "):
            current = line[4:].strip()
        elif current and line.strip().startswith(label):
            found[current] = line.split(label, 1)[1].strip()
            current = None
    return found


def parse_masked_candidates(report: str) -> list[str]:
    """아형 마스킹 리포트에서 후보 유전자 목록을 뽑아낸다."""
    for line in report.splitlines():
        if line.startswith("SUBTYPE_MASKED_CANDIDATES:"):
            return [g.strip() for g in line.split(":", 1)[1].split(",") if g.strip()]
    return []


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
    parser.add_argument("--skip-organoid", action="store_true", help="오가노이드 관점 생략 (세포주 관점만)")
    parser.add_argument(
        "--csv",
        default=None,
        help="후보 표를 CSV로도 저장할 경로 (기본: 저장하지 않고 화면에만 출력)",
    )
    args = parser.parse_args()

    print("=" * 78)
    print(f"사용자 질문: {args.cancer_type}에서 {args.mutation} 변이 기반 합성치사 후보를 찾아줘")
    print("=" * 78)
    print("")

    # ---- 1단계: 실제 DepMap 데이터로 후보 발굴 -------------------------------------------
    # 후보표는 화면에 출력하는 것이 기본이다. --csv를 준 경우에만 파일로 남기고, 그렇지 않으면
    # 임시 파일에 받아서 읽은 뒤 지운다 (작업 디렉터리에 산출물을 남기지 않기 위해).
    import tempfile

    import pandas as pd

    table_path = args.csv or os.path.join(tempfile.mkdtemp(prefix="biomni_sl_"), "candidates.csv")
    report = discover_synthetic_lethal_candidates(
        cancer_type=args.cancer_type,
        target_mutation=args.mutation,
        top_n=args.top_n,
        output_csv_path=table_path,
    )
    print(report)

    if report.startswith("FAILURE") or not os.path.exists(table_path):
        print("\n후보 발굴에 실패했습니다. 위의 실패 사유를 확인하세요.")
        return

    table = pd.read_csv(table_path)
    if not args.csv:
        os.remove(table_path)
        os.rmdir(os.path.dirname(table_path))
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
    print_interpretation_guide(table, targets)

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
        print(f"  -> 남은 실험 대상: {', '.join(clean) if clean else '없음'}")
        print("")

    if not clean:
        print("모든 후보가 교란요인으로 설명됩니다. 실험을 시작하기 전에 해당 바이오마커로")
        print("층화한 뒤 discover_synthetic_lethal_candidates를 다시 실행하세요.")
        return

    # ---- 5단계: 오가노이드 관점 -------------------------------------------------------------
    masked_genes: list[str] = []
    transfer: dict = {}
    if not args.skip_organoid:
        print("=" * 78)
        print("오가노이드 관점 (1) — 세포주 통합 분석이 아형 평균화로 놓친 후보")
        print("=" * 78)
        masked_report = discover_subtype_masked_sl_candidates(
            cancer_type=args.cancer_type,
            target_mutation=args.mutation,
            top_n=args.top_n,
        )
        print(masked_report)
        masked_genes = parse_masked_candidates(masked_report)

        print("")
        print("=" * 78)
        print("오가노이드 관점 (2) — 각 후보가 오가노이드에서도 검출될 것인가")
        print("=" * 78)
        assessed = list(dict.fromkeys(clean + masked_genes))
        transfer_report = assess_organoid_transferability(
            candidate_genes=assessed,
            target_mutation=args.mutation,
            cancer_type=args.cancer_type,
        )
        print(transfer_report)
        transfer = parse_labelled_block(transfer_report, "VERDICT:")

        # ---- 통합 요약 ---------------------------------------------------------------------
        print("")
        print("=" * 78)
        print("통합 요약 — 세포주 관점 vs 오가노이드 관점")
        print("=" * 78)
        print(f"{'유전자':<12}{'발굴 관점':<22}{'교란 판정':<18}오가노이드 전이성")
        print("-" * 78)
        for gene in assessed:
            origin = "세포주 (통합 검정)" if gene in clean else "오가노이드 (아형 마스킹)"
            confound = verdicts.get(gene, "-").split(" - ")[0]
            organoid = transfer.get(gene, "-").split(" - ")[0]
            print(f"{gene:<12}{origin:<22}{confound:<18}{organoid}")
        print("")
        print("  · 세포주 관점 후보: 근거는 가장 강하지만 2D 배양·정상 대조군 부재의 한계를 안습니다.")
        print("  · 오가노이드 관점 후보: 근거는 약하지만, 세포주 패널로는 원리적으로 판정할 수 없어")
        print("    오가노이드에서만 답이 나옵니다.")
        print("  · ORGANOID-MASKED 후보는 표준 배지에서 위음성이 납니다. 배지 수정 없이 실험하지 마세요.")
        print("")

    # ---- 6단계: 최우선 후보의 실험 프로토콜 ------------------------------------------------
    top = targets[targets["gene"] == clean[0]].iloc[0]
    print_experiment_plan(top["gene"], top, args.mutation)
    print(f"교란요인 판정: {verdicts.get(top['gene'], 'N/A')}")

    if not args.skip_organoid:
        print("")
        print(
            design_organoid_sl_experiment(
                candidate_gene=top["gene"],
                target_mutation=args.mutation,
                cancer_type=args.cancer_type,
            )
        )
        print("")
        print(f"오가노이드 전이성 판정: {transfer.get(top['gene'], 'N/A')}")

    if args.csv:
        print(f"\n전체 후보 표를 {args.csv}에 저장했습니다.")


if __name__ == "__main__":
    main()
