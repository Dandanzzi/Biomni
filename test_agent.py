"""유방암 BRCA1/BRCA2 합성치사 후보 발굴 — 다채널 통계 방식.

Noh et al., New Biotechnology 94:184-191 (2026) 리뷰가 지적하는 문제를 반영해 다시 짰습니다.

이전 방식(discover_synthetic_lethal_candidates)은 단일 채널 + genome-wide FDR 게이트였고,
유방암 BRCA에서 후보를 0개 냈습니다. 교과서적인 BRCA-PARP1조차 q=0.998로 탈락했습니다.
리뷰 3장이 말하듯 합성치사 발굴은 극단적 class imbalance 아래의 '순위' 문제이지 임계값
문제가 아니며, 임계-후-집계 방식은 성능을 과대평가하거나(모델의 경우) 아무것도 못 찾습니다.

바뀐 것:

  1. 단일 채널 -> DAISY식 3채널 (리뷰 2.1장, Jerby-Arnon Cell 2014)
       essentiality    driver 결손주에서 더 필수적인가
       co-expression   driver와 발현이 함께 움직이는가 (기능적 근접성)
       co-inactivation 둘 다 잃은 세포주가 드문가 (survival of the fittest)
     각 채널은 개별로 약하지만 오차가 서로 독립이라, 순위를 합치면 신호가 남습니다.

  2. FDR 게이트 -> 순위 + 알려진 상호작용 벤치마크 (리뷰 3장)
     q<0.25로 자르는 대신 전체를 순위로 내고, 이미 확립된 파트너(POLQ, PARP1/2)가
     몇 위에 오는지를 함께 출력합니다. 그 순위가 이 실행을 믿을지 말지의 근거입니다.

  3. essentiality floor를 gene effect 컷 -> 의존 '확률' (리뷰 5장)
     Project Score의 Bayes factor와 같은 취지. 상시 필수유전자가 목록을 덮는 것을 막습니다.

  4. druggability를 사후 필터 -> 점수 단계에 포함 (리뷰 5장)
     Broad Repurposing Hub의 임상 단계를 후보와 함께 답니다.

  5. 유방암 단독 -> pan-cancer 기본값
     BRCA 결손 표현형은 조직 비특이적인데 유방암 패널만 쓰면 결손주가 6종뿐입니다.
     pan-cancer로 50종 이상이 모여야 어느 채널이든 검정력이 생깁니다.

    python test_agent.py                          # BRCA1+BRCA2 / pan-cancer (기본값)
    python test_agent.py --cancer-type "Breast Cancer"   # 유방암으로 좁혀서 (검정력 낮음)
    python test_agent.py --mutation BRCA1         # 단일 driver
    python test_agent.py --skip-pubmed            # 문헌 검증 생략
"""

import argparse

from biomni.tool.sl_multichannel import discover_sl_multichannel
from biomni.tool.synthetic_lethality import (
    check_dependency_confounders,
    validate_sl_candidates_with_pubmed,
)


def parse_candidate_line(report: str, label: str) -> list[str]:
    """'LABEL: GENE1, GENE2, ...' 한 줄에서 유전자 목록을 뽑아낸다.

    각 발굴 도구는 사람이 읽는 리포트 끝에 이런 줄을 하나 남긴다. 표를 다시 파싱하지 않고
    이 줄만 읽으면 되므로, 도구의 출력 서식이 바뀌어도 파이프라인이 깨지지 않는다.
    """
    prefix = f"{label}:"
    for line in report.splitlines():
        if line.startswith(prefix):
            genes = [g.strip() for g in line.split(":", 1)[1].split(",") if g.strip()]
            return [g for g in genes if g != "(none)"]
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


def describe_cell_line_panel(mutation: str, cancer_type: str, max_shown: int = 6):
    """통계를 낸 그 패널에서 변이/야생형 세포주를 그대로 읽어온다.

    run_agent.py는 췌장암 세포주 이름(MIA PaCa-2, PANC-1, BxPC-3)을 코드에 박아 두었다.
    질문만 유방암으로 바꾸면 그 목록이 통째로 거짓이 되므로, 여기서는 DepMap에서 직접
    읽는다. 화면에 뜬 세포주가 곧 근거가 나온 세포주라 목록과 통계가 어긋날 수 없다.
    """
    try:
        from biomni.tool.synthetic_lethality import (
            _annotate_mutation_status,
            _load_depmap,
            _select_cancer_models,
        )

        bundle = _load_depmap(None)
        cohort, _ = _select_cancer_models(bundle["model"], cancer_type)
        cohort = cohort[cohort["ModelID"].isin(bundle["gene_effect"].index)]
        cohort, source = _annotate_mutation_status(cohort, mutation, bundle["data_lake_path"], None)
    except Exception:  # 패널 표시는 부가 정보라, 실패해도 계획 출력은 계속한다
        return None

    mutant = cohort[cohort["MutationStatus"] == "MUT"]
    described = [
        f"{row['StrippedCellLineName']} ({row['Variant']})" if row["Variant"] else row["StrippedCellLineName"]
        for _, row in mutant.iterrows()
    ]
    wildtype = cohort.loc[cohort["MutationStatus"] == "WT", "StrippedCellLineName"].tolist()
    return described[:max_shown], wildtype[:max_shown], len(mutant), len(wildtype), source


def print_experiment_plan(gene: str, mutation: str, cancer_type: str) -> None:
    print("=" * 78)
    print(f"실험 계획: {mutation} 변이 의존적 {gene} 합성치사 검증")
    print("=" * 78)
    print(f"  가설: {mutation} 변이 {cancer_type} 세포는 {gene} 결손에 선택적으로 취약하다.")
    print("")
    print("  1) 세포주 패널 (DepMap에서 실제로 스크리닝된 계열)")
    panel = describe_cell_line_panel(mutation, cancer_type)
    if panel is None:
        print("     - DepMap 패널을 읽지 못했습니다. 위 통계를 낸 세포주 목록을 직접 확인하세요.")
    else:
        mutant_lines, wildtype_lines, n_mutant, n_wildtype, source = panel
        print(f"     - {mutation} 변이군 ({n_mutant}종): {', '.join(mutant_lines)}")
        print(f"     - {mutation} 야생형군 ({n_wildtype}종): {', '.join(wildtype_lines)}")
        print(f"     * 통계를 낸 그 패널에서 그대로 읽은 목록입니다 (출처: {source}).")
        print("     * 어느 한쪽 군이 10종 미만이면 그 군만으로 대조하지 말고 isogenic 계통을 함께 쓸 것.")
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
    print("     - 비형질전환 정상 세포에서 동등한 의존성이 나오면 therapeutic window 없음 -> No-go")
    print("")


def main() -> None:
    parser = argparse.ArgumentParser(description="다채널 통계로 합성치사 후보를 발굴한다")
    # 여기가 프롬프트다.
    parser.add_argument("--cancer-type", default="pan-cancer")
    parser.add_argument("--mutation", default="BRCA1,BRCA2", help="driver. 쉼표로 여러 개를 주면 결손을 풀링한다")
    parser.add_argument(
        "--known-sl",
        default="PARP1,PARP2,POLQ,RAD51,ATR",
        help="이미 확립된 파트너. 순위 보정 확인용이며 발굴 대상이 아니다",
    )
    parser.add_argument("--top-k", type=int, default=20, help="보고할 후보 수")
    parser.add_argument("--require-channels", type=int, default=3, help="몇 개 채널에서 점수가 나야 하는가")
    parser.add_argument("--no-lof-only", dest="lof_only", action="store_false", help="missense VUS도 결손으로 취급")
    parser.add_argument("--skip-pubmed", action="store_true", help="문헌 검증 생략 (오프라인 실행)")
    parser.add_argument("--validate-top", type=int, default=6, help="문헌·교란 검사로 넘길 상위 후보 수")
    args = parser.parse_args()

    genes = [g.strip().upper() for g in args.mutation.split(",") if g.strip()]
    known = [g.strip().upper() for g in args.known_sl.split(",") if g.strip()]

    print("=" * 78)
    print(f"사용자 질문: {args.cancer_type}에서 {'/'.join(genes)} 결손 기반 합성치사 후보를 찾아줘")
    print("=" * 78)
    print("")

    # ---- 1단계: 다채널 발굴 -------------------------------------------------------------------
    report = discover_sl_multichannel(
        driver_genes=genes,
        cancer_type=args.cancer_type,
        known_sl_partners=known,
        lof_only=args.lof_only,
        top_k=args.top_k,
        require_channels=args.require_channels,
    )
    print(report)
    if report.startswith("FAILURE"):
        return
    candidates = parse_candidate_line(report, "MULTICHANNEL_CANDIDATES")
    if not candidates:
        print("\n순위에 오른 후보가 없습니다. --require-channels 를 낮춰 보세요.")
        return

    validate = candidates[: args.validate_top]
    print("")
    print(f"=> 문헌·교란 검사로 넘길 상위 후보: {', '.join(validate)}")
    print("   (순위는 통계 신호이고, 아래 두 단계가 그 신호를 반증하는 자리입니다)")
    print("")

    # ---- 2단계: 문헌 검증 ---------------------------------------------------------------------
    if not args.skip_pubmed:
        print(
            validate_sl_candidates_with_pubmed(
                disease=args.cancer_type if args.cancer_type != "pan-cancer" else "cancer",
                mutated_gene=genes[0],
                candidate_genes=validate,
                max_papers_per_gene=6,
                include_abstracts=False,
            )
        )
        print("")

    # ---- 3단계: 교란요인 반증 검사 -------------------------------------------------------------
    confounder_report = check_dependency_confounders(
        target_mutation=genes[0],
        candidate_genes=validate,
        cancer_type=args.cancer_type,
    )
    print(confounder_report)
    verdicts = parse_confounder_verdicts(confounder_report)

    clean = [g for g in validate if not verdicts.get(g, "").startswith("CONFOUNDED")]
    dropped = [g for g in validate if g not in clean]
    if dropped:
        print("")
        print("=" * 78)
        print("교란요인으로 제외된 후보")
        print("=" * 78)
        for gene in dropped:
            print(f"  {gene}: {verdicts[gene]}")
        print("")
    if not clean:
        print("남은 후보가 없습니다. 해당 바이오마커로 층화한 뒤 다시 실행하세요.")
        return

    # ---- 4단계: 최우선 후보의 실험 계획 --------------------------------------------------------
    print_experiment_plan(clean[0], genes[0], args.cancer_type)
    print(f"교란요인 판정: {verdicts.get(clean[0], 'N/A')}")
    print("")
    print("=" * 78)
    print(f"최종 — {'/'.join(genes)} 결손 합성치사 상위 {min(len(clean), 5)}순위")
    print("=" * 78)
    for index, gene in enumerate(clean[:5], start=1):
        print(f"  {index}. {gene:<10}{verdicts.get(gene, '-').split(' - ')[0]}")
    print("")
    print("  이 순위는 통계 채널의 합의일 뿐 검증된 상호작용이 아닙니다. 위 STEP 5의")
    print("  보정 결과(알려진 파트너가 몇 위였는지)를 먼저 보고 신뢰도를 판단하세요.")


if __name__ == "__main__":
    main()
