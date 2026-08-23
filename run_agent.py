"""췌장암 KRAS 합성치사: 세포주 관점과 오가노이드 관점을 나란히 산출한다.

LLM을 쓰지 않는 결정론적 파이프라인이라 같은 입력이면 항상 같은 결과가 나옵니다.
(에이전트가 스스로 도구를 계획하는 버전은 run_sl_agent.py --mode agent)

두 관점을 모두 출력합니다.

  [세포주 관점]  DepMap CRISPR 47개 췌장암 세포주에서 통계적으로 발굴한 후보.
                 근거는 가장 강하지만, 2D 배양·아형 평균화·유전형 라벨 뭉뚱그리기·정상 대조군
                 부재라는 한계를 안습니다.
  [오가노이드 관점] 세포주 패널이 원리적으로 답할 수 없는 후보를 세 축으로 나눠 발굴합니다.
                 (1) 아형 평균화로 상쇄된 것   (2) KRAS allele을 하나로 묶어 희석된 것
                 (3) 혈청이 지질을 공급해 2D에서는 보이지 않는 것
                 그리고 각 후보가 오가노이드 배지에서 가려질지(위음성) 여부까지 판정합니다.

(2)와 (3)은 tumour-derived organoid biobank 논문(Nature 2026; 환자 유래 오가노이드 256개,
그중 162개 genome-wide CRISPR)의 두 발견을 로컬 DepMap에서 검증 가능한 형태로 옮긴 것입니다.
그 논문은 세포주 DepMap과 공유하는 core fitness 654개 외에 오가노이드 특이적 core fitness
97개를 보고했고(steroid/cholesterol/isoprenoid 생합성에 집중), KRAS 의존성이 allele에 따라
갈린다는 것을 보였습니다(G12 계열은 KRAS·EGFR·PTPN11 의존, Q61H는 EGFR 억제에 무반응).

    python run_agent.py                      # 췌장암 / KRAS (기본값)
    python run_agent.py --mutation TP53      # 다른 driver 변이
    python run_agent.py --skip-pubmed        # 문헌 검증 생략
    python run_agent.py --skip-organoid      # 세포주 관점만
    python run_agent.py --skip-serum-axis    # 지질/혈청 축만 생략 (가장 오래 걸리는 단계)
    python run_agent.py --focus-gene SREBF1  # 상세 실험 계획을 세울 후보를 직접 지정
"""

import argparse
import os

from biomni.tool.organoid_sl import (
    assess_organoid_transferability,
    design_organoid_sl_experiment,
    discover_allele_resolved_sl_candidates,
    discover_serum_masked_dependencies,
    discover_subtype_masked_sl_candidates,
    rank_organoid_sl_candidates,
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
    print("  선택도             변이의존% - 전체의존%. genotype 특이성의 핵심 지표입니다.")
    print("  * 실제 우선순위는 선택도 하나가 아니라 효과크기·치료 window·FDR을 합친 점수로")
    print("    정합니다 (select_experiment_targets). 성분은 아래 선별 표에 실려 있습니다.")
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


def add_priority_columns(frame):
    """우선순위 점수와 그 성분들을 컬럼으로 붙인다 (원본은 건드리지 않는다).

    선별 통과 후보든, --focus-gene으로 나중에 끌어올린 후보든 같은 함수를 거치게 해서
    두 경로가 서로 다른 점수 체계를 갖는 일이 없도록 한다.
    """
    scored = frame.copy()
    # 선택도(변이 세포주 의존 비율 - 전체 세포주 의존 비율)가 클수록 genotype 특이적이다.
    scored["selectivity_gap"] = scored["pct_mutant_dependent"] - scored["pct_all_lines_dependent"]
    # 효과크기: 변이군과 야생형군의 gene effect 차이. 음수여야 하므로 절대값을 쓴다.
    scored["effect_component"] = scored["effect_difference"].abs().clip(upper=0.5) * 40
    # 치료 window: 전체 세포주 의존 비율이 20%를 넘는 만큼 깎는다. 20% 아래는 벌점 없음.
    scored["window_penalty"] = (scored["pct_all_lines_dependent"] - 20.0).clip(lower=0) * 0.5
    # 다중검정: q가 클수록 깎는다. q=0.05면 -5점, 게이트 상한 q=0.10이면 -10점.
    scored["fdr_penalty"] = scored["q_value"].clip(upper=MAX_Q_VALUE) * 100
    # 야생형 세포주도 의존하면 genotype 선택적이지 않다.
    scored["wildtype_penalty"] = scored["pct_wildtype_dependent"] * 0.5
    scored["priority_score"] = (
        scored["selectivity_gap"]
        + scored["effect_component"]
        - scored["window_penalty"]
        - scored["fdr_penalty"]
        - scored["wildtype_penalty"]
    )
    return scored


def select_experiment_targets(table, top_k: int = 5):
    """통계 결과표에서 wet-lab 검증 우선순위를 규칙 기반으로 선별한다.

    선별 자체는 세 개의 하드 게이트(pan-essential / FDR / 변이주 의존 비율)이고, 통과한
    후보의 '순서'는 아래 네 성분의 합으로 정한다. 선택도 하나로 정렬하면 pan-essential
    쪽으로 새는 후보와 FDR에 간신히 걸린 후보가 위로 올라오는데, 두 결함 모두 선택도에는
    나타나지 않기 때문이다.

    각 성분을 컬럼으로 남기므로 순위에 이의가 있으면 어느 성분 때문인지 표에서 바로
    짚을 수 있다. 파라로그 교란은 여기서 계산하지 않는다 - 발현 상관이 필요하고,
    check_dependency_confounders와 rank_organoid_sl_candidates가 이미 그 검사를 맡는다.
    이 단계는 문헌도 오가노이드 축도 보지 못한다는 점을 전제로 읽어야 한다.
    """
    eligible = table[
        (table["pct_all_lines_dependent"] < MAX_PAN_ESSENTIAL_PCT)
        & (table["q_value"] < MAX_Q_VALUE)
        & (table["pct_mutant_dependent"] >= MIN_MUTANT_DEPENDENT_PCT)
    ]
    eligible = add_priority_columns(eligible)
    eligible = eligible.sort_values(["priority_score", "effect_difference"], ascending=[False, True])
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


def parse_literature_scores(report: str) -> dict:
    """PubMed 리포트에서 유전자별 0-100 문헌 점수를 뽑아낸다."""
    scores = {}
    current = None
    for line in report.splitlines():
        if line.startswith("### "):
            current = line[4:].strip()
        elif current and "Literature support score:" in line:
            scores[current] = int(line.split("Literature support score:", 1)[1].split("/")[0].strip())
            current = None
    return scores


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


def parse_ranking_rows(report: str) -> list[dict]:
    """rank_organoid_sl_candidates의 ORGANOID_SL_ROW 줄을 순위 순서대로 읽는다."""
    rows = []
    for line in report.splitlines():
        if not line.startswith("ORGANOID_SL_ROW:"):
            continue
        parts = line.split(":", 1)[1].strip().split("|")
        if len(parts) != 7:
            continue
        gene, score, klass, reading, difference, pct_all, normal_tpm = parts
        rows.append(
            {
                "gene": gene,
                "score": float(score),
                "class": klass,
                "reading": reading,
                "difference": float(difference) if difference else None,
                "pct_all": float(pct_all) if pct_all else None,
                "normal_tpm": float(normal_tpm) if normal_tpm else None,
            }
        )
    return rows


def print_final_ranking(ranking_rows, targets, verdicts, mutation: str, top_k: int = 5) -> None:
    """두 관점의 상위 후보를 나란히 출력한다.

    같은 데이터를 두 가지 방식으로 정렬한 결과이므로, 1순위가 갈리는 것이 정상이다.
    갈리는 이유를 마지막에 한 줄로 적어 두 표를 어떻게 읽어야 하는지 남긴다.
    """
    print("")
    print("=" * 78)
    print(f"최종 종합 — 두 관점의 상위 {top_k}순위")
    print("=" * 78)

    print("")
    print(f"[오가노이드 관점] 합성치사에 가장 가까운 후보 {top_k}개")
    print("  기준: rank_organoid_sl_candidates. 유전형 선택성 게이트 통과 + 오가노이드가")
    print("        실제로 답을 더 잘 주는가(부착/니치)를 함께 봅니다.")
    print("")
    organoid_top = []
    if ranking_rows:
        organoid_top = [r for r in ranking_rows if r["class"] == "SL-CANDIDATE"][:top_k]
    if not organoid_top:
        print("  (오가노이드 순위를 산출하지 않았습니다 — --skip-organoid 이거나 통과 후보 없음)")
    else:
        print(f"  {'순위':<5}{'유전자':<12}{'점수':>7}{'MUT-WT':>9}{'%전체의존':>10}  {'판정':<16}오가노이드 판독")
        print("  " + "-" * 74)
        for index, row in enumerate(organoid_top, start=1):
            difference = f"{row['difference']:+.3f}" if row["difference"] is not None else "-"
            pct_all = f"{row['pct_all']:.0f}%" if row["pct_all"] is not None else "-"
            print(
                f"  {index:<5}{row['gene']:<12}{row['score']:>7.1f}{difference:>9}{pct_all:>10}  "
                f"{row['class']:<16}{row['reading']}"
            )

    print("")
    print(f"[세포주 관점] 합성치사에 가장 가까운 후보 {top_k}개")
    print("  기준: select_experiment_targets의 우선순위 점수 정렬,")
    print("        여기에 check_dependency_confounders의 반증 판정을 덧붙였습니다.")
    print("")
    # 교란요인으로 기각된 후보는 순위에서 뺀다. 오가노이드 표가 SL-CANDIDATE만 싣는 것과
    # 같은 이유이고, 무엇보다 파이프라인이 스스로 기각한 후보를 '합성치사에 가장 가까운'
    # 표의 1위로 올리면 그 앞 단계의 반증 검사를 무효로 만든다.
    ranked = targets.sort_values("priority_score", ascending=False)
    confounded = [g for g in ranked["gene"] if verdicts.get(g, "").startswith("CONFOUNDED")]
    cellline_top = ranked[~ranked["gene"].isin(confounded)].head(top_k)
    print(
        f"  {'순위':<5}{'유전자':<12}{'점수':>7}{'선택도':>8}{'효과':>7}"
        f"{'window':>8}{'FDR':>6}  교란 판정"
    )
    print("  " + "-" * 74)
    for index, (_, row) in enumerate(cellline_top.iterrows(), start=1):
        verdict = verdicts.get(row["gene"], "-").split(" - ")[0]
        print(
            f"  {index:<5}{row['gene']:<12}{row['priority_score']:>7.1f}{row['selectivity_gap']:>8.1f}"
            f"{row['effect_component']:>7.1f}{-row['window_penalty']:>8.1f}{-row['fdr_penalty']:>6.1f}  {verdict}"
        )
    if confounded:
        print(f"  (교란요인으로 제외: {', '.join(confounded)})")
    if len(cellline_top) < top_k:
        print(f"  (선별 기준을 통과하고 교란 검사까지 살아남은 후보가 {len(cellline_top)}개뿐입니다)")
    print("  · 점수 = 선택도 + 효과 - window - FDR - 야생형의존. 성분을 그대로 실어 두었으므로")
    print("    순위에 이의가 있으면 어느 성분 때문인지 짚을 수 있습니다.")
    print("  · 다만 이 단계는 파라로그 교란도 문헌도 오가노이드 축도 보지 못합니다.")
    print("    그래서 두 표의 순위가 갈립니다.")

    print("")
    if organoid_top and len(cellline_top):
        organoid_first = organoid_top[0]["gene"]
        cellline_first = cellline_top.iloc[0]["gene"]
        if organoid_first == cellline_first:
            print(f"  두 관점의 1순위가 {organoid_first}로 일치합니다.")
        else:
            print(f"  두 관점의 1순위가 다릅니다: 오가노이드 {organoid_first} vs 세포주 {cellline_first}.")
            print("  아래 표는 DepMap 통계만으로 매긴 순위입니다. 파라로그 교란·문헌·오가노이드가")
            print("  실제로 무엇을 더 측정해 주는지는 그 시점에 아직 계산되지 않았습니다.")
            print(f"  오가노이드에 예산을 쓸 것인지 묻는 질문이라면 위쪽 표({organoid_first})가 답입니다.")
    print("")


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
        "--skip-serum-axis", action="store_true", help="지질/혈청 축 생략 (permutation 검정이 오래 걸림)"
    )
    parser.add_argument(
        "--serum-permutations",
        type=int,
        default=2000,
        help="지질 프로그램 농축 검정의 permutation 횟수 (기본 2000)",
    )
    parser.add_argument(
        "--focus-gene",
        default=None,
        help="상세 실험 계획을 세울 후보 유전자 (기본: 선별 1순위 후보. 예: --focus-gene SREBF1)",
    )
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

    # --focus-gene를 준 경우, 선별 상위권 밖이더라도 후보표에 있으면 분석 대상에 끌어올린다.
    # (문헌 검증·교란 검사·오가노이드 전이성 판정까지 모두 그 유전자를 포함해서 돌아야 한다.)
    focus_gene = args.focus_gene.strip().upper() if args.focus_gene else None
    if focus_gene and focus_gene not in set(targets["gene"]):
        extra = table[table["gene"] == focus_gene].copy()
        if extra.empty:
            print("")
            print(f"  --focus-gene {focus_gene}: 상위 {args.top_n}개 후보표에 없습니다.")
            print("  --top-n을 늘려 다시 실행하거나, 위 표에 있는 유전자를 지정하세요.")
            return
        targets = pd.concat([targets, add_priority_columns(extra)], ignore_index=True)
        print(f"  --focus-gene {focus_gene}: 선별 상위권 밖이지만 지정에 따라 분석 대상에 추가합니다.")
    print("")
    print_interpretation_guide(table, targets)

    print("")
    print(
        f"{'순위':<6}{'유전자':<12}{'점수':>7}{'변이군':>9}{'야생형':>9}{'차이':>9}{'q':>9}"
        f"{'변이의존':>9}{'전체의존':>9}"
    )
    print("-" * 78)
    for index, (_, row) in enumerate(targets.iterrows(), start=1):
        print(
            f"{index:<6}{row['gene']:<12}{row['priority_score']:>7.1f}"
            f"{row['mutant_mean_effect']:>9.3f}{row['wildtype_mean_effect']:>9.3f}"
            f"{row['effect_difference']:>9.3f}{row['q_value']:>9.4f}"
            f"{row['pct_mutant_dependent']:>8.0f}%{row['pct_all_lines_dependent']:>8.0f}%"
        )
    print("  점수 = 선택도 + 효과크기 - window벌점 - FDR벌점 - 야생형의존벌점 (select_experiment_targets)")

    candidate_genes = targets["gene"].tolist()
    print("")
    print(f"=> 문헌 검증 및 실험 대상: {', '.join(candidate_genes)}")
    print("   (이전 버전의 'Gene_B'와 달리 모두 실존하는 HUGO 심볼이며 그대로 sgRNA 주문에 사용할 수 있습니다)")
    print("")

    # ---- 3단계: PubMed 검증 ----------------------------------------------------------------
    literature_scores: dict[str, int] = {}
    if not args.skip_pubmed:
        pubmed_report = validate_sl_candidates_with_pubmed(
            disease=args.cancer_type,
            mutated_gene=args.mutation,
            candidate_genes=candidate_genes,
            max_papers_per_gene=6,
            include_abstracts=False,
        )
        print(pubmed_report)
        literature_scores = parse_literature_scores(pubmed_report)
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
    # 세 개의 독립적인 축으로 나눠서 발굴한다. 각 축은 세포주 패널이 "아직 못 찾은" 것이 아니라
    # "구조적으로 찾을 수 없는" 것을 겨냥한다 — 그래서 오가노이드가 확인 단계가 아니라 유일한
    # 실험이 되는 지점이다.
    origins: dict[str, list[str]] = {}
    transfer: dict = {}
    organoid_pick: str | None = None
    ranking_rows: list[dict] = []
    if not args.skip_organoid:
        print("=" * 78)
        print("오가노이드 관점 (1) — 아형 평균화로 상쇄된 후보")
        print("=" * 78)
        masked_report = discover_subtype_masked_sl_candidates(
            cancer_type=args.cancer_type,
            target_mutation=args.mutation,
            top_n=args.top_n,
        )
        print(masked_report)
        for gene in parse_candidate_line(masked_report, "SUBTYPE_MASKED_CANDIDATES"):
            origins.setdefault(gene, []).append("아형")

        print("")
        print("=" * 78)
        print(f"오가노이드 관점 (2) — '{args.mutation} 변이' 라벨이 allele을 뭉뚱그려 희석한 후보")
        print("=" * 78)
        allele_report = discover_allele_resolved_sl_candidates(
            cancer_type=args.cancer_type,
            target_mutation=args.mutation,
            top_n=args.top_n,
        )
        print(allele_report)
        for gene in parse_candidate_line(allele_report, "ALLELE_RESOLVED_CANDIDATES"):
            origins.setdefault(gene, []).append("allele")

        if not args.skip_serum_axis:
            print("")
            print("=" * 78)
            print("오가노이드 관점 (3) — 혈청이 지질을 공급해 2D에서는 보이지 않는 후보")
            print("=" * 78)
            serum_report = discover_serum_masked_dependencies(
                cancer_type=args.cancer_type,
                target_mutation=args.mutation,
                n_permutations=args.serum_permutations,
                top_n=args.top_n,
            )
            print(serum_report)
            for gene in parse_candidate_line(serum_report, "SERUM_MASKED_CANDIDATES"):
                origins.setdefault(gene, []).append("지질")

        print("")
        print("=" * 78)
        print("오가노이드 관점 (4) — 각 후보가 오가노이드 배지에서도 검출될 것인가")
        print("=" * 78)
        assessed = list(dict.fromkeys(clean + list(origins)))
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
        print(f"{'유전자':<12}{'발굴 축':<20}{'교란 판정':<18}오가노이드 전이성")
        print("-" * 78)
        for gene in assessed:
            axes = origins.get(gene, [])
            if gene in clean:
                axes = ["세포주"] + axes
            confound = verdicts.get(gene, "-").split(" - ")[0]
            organoid = transfer.get(gene, "-").split(" - ")[0]
            print(f"{gene:<12}{'+'.join(axes):<20}{confound:<18}{organoid}")
        print("")

        # ---- 오가노이드 관점 (5) — 후보 간 순위와 결론 ---------------------------------------
        # 여기까지의 도구들은 유전자별 판정만 내리고 후보끼리 비교하지 않는다. 어느 후보에
        # 오가노이드를 쓸 것인가라는 결론은 지금까지 사람이 리포트를 읽고 내려야 했는데,
        # rank_organoid_sl_candidates가 그 판단을 명시적 규칙으로 대신한다. 규칙과 임계값은
        # 도구가 출력에 그대로 찍으므로, 결론에 이의가 있으면 임계값을 짚어 반박할 수 있다.
        print("=" * 78)
        print("오가노이드 관점 (5) — 후보 간 순위와 결론 (규칙 기반)")
        print("=" * 78)
        axis_labels = {"아형": "subtype", "allele": "allele", "지질": "lipid"}
        axis_origins = {}
        for gene in assessed:
            axes = ["cell-line"] if gene in clean else []
            axes += [axis_labels.get(a, a) for a in origins.get(gene, [])]
            axis_origins[gene] = axes
        ranking_report = rank_organoid_sl_candidates(
            candidate_genes=assessed,
            target_mutation=args.mutation,
            cancer_type=args.cancer_type,
            axis_origins=axis_origins,
            literature_scores=literature_scores or None,
            q_values={row["gene"]: row["q_value"] for _, row in table.iterrows()},
        )
        print(ranking_report)
        ranked_top = parse_candidate_line(ranking_report, "ORGANOID_TOP_PICK")
        organoid_pick = ranked_top[0] if ranked_top else None
        ranking_rows = parse_ranking_rows(ranking_report)
        print("")

    # ---- 6단계: 최우선 후보의 실험 프로토콜 ------------------------------------------------
    # 상세 계획 대상은 세 단계로 정한다.
    #   1) --focus-gene을 준 경우 그 유전자 (사용자 지정이 최우선)
    #   2) 오가노이드 순위 도구가 고른 ORGANOID_TOP_PICK (규칙 기반 결론)
    #   3) 둘 다 없으면 교란 검사를 통과한 선별 1순위
    # 2)를 3)보다 앞에 두는 이유: 선별 1순위는 선택도 하나로 정렬한 결과라 유전형 선택성·
    # 치료 window·파라로그 교란을 반영하지 않는다.
    detail_gene = focus_gene or organoid_pick or clean[0]
    if not focus_gene and organoid_pick and organoid_pick != clean[0]:
        print(f"※ 상세 계획 대상을 선별 1순위({clean[0]})가 아니라 오가노이드 순위 도구의 결론")
        print(f"   ORGANOID_TOP_PICK={organoid_pick}으로 잡습니다. 근거는 바로 위 랭킹 리포트에 있습니다.")
        print("")
    if detail_gene not in set(targets["gene"]):
        extra_row = table[table["gene"] == detail_gene]
        if extra_row.empty:
            print(f"{detail_gene}의 통계 행이 후보표에 없어 상세 계획을 건너뜁니다.")
            return
        targets = pd.concat([targets, add_priority_columns(extra_row)], ignore_index=True)
    if focus_gene and focus_gene not in clean:
        print("=" * 78)
        print(f"주의: {focus_gene}는 교란요인 검사를 통과하지 못했습니다 ({verdicts.get(focus_gene, 'N/A')}).")
        print("지정에 따라 실험 계획은 출력하지만, 위 판정을 먼저 해소한 뒤 진행하세요.")
        print("=" * 78)
        print("")
    top = targets[targets["gene"] == detail_gene].iloc[0]
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

    print_final_ranking(ranking_rows, targets, verdicts, args.mutation)

    if args.csv:
        print(f"\n전체 후보 표를 {args.csv}에 저장했습니다.")


if __name__ == "__main__":
    main()
