# Biomni Cancer Vulnerability Toolkit

**합성치사 후보 발굴 → 반증 검사 → 오가노이드 검증 → 예측 채점까지, Biomni 호환 도구 생태계**

DepMap CRISPR 의존성으로 후보를 발굴하고, PubMed·STRING으로 교차 검증하고, **"이 후보가 정말 driver 변이 때문인가"를
능동적으로 반증**한 뒤, 환자유래 오가노이드(PDO) 실측 데이터를 읽어 들여 **예측이 맞았는지 채점**합니다.

> 이 도구는 "점수를 잘 내는 예측기"가 아니라 **"어떤 후보를 왜 먼저 검증해야 하는가"에 답하는 의사결정 계층**입니다.
> 그래서 지지 근거만 모으지 않고 pan-essentiality·파라로그 손실·lineage 교란·문헌상 반증을 적극적으로 찾아 후보를 탈락시키며,
> 마지막 단계는 **"예측이 틀렸다"고 말할 수 있도록** 설계되어 있습니다.

---

## 기술 스택

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![pandas](https://img.shields.io/badge/pandas-2.x-150458?logo=pandas&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-2.x-013243?logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-1.15-8CAAE6?logo=scipy&logoColor=white)
![matplotlib](https://img.shields.io/badge/matplotlib-3.10-11557C?logo=python&logoColor=white)
![seaborn](https://img.shields.io/badge/seaborn-0.13-4C72B0)
![Ruff](https://img.shields.io/badge/lint-Ruff-D7FF64?logo=ruff&logoColor=black)

![Biomni](https://img.shields.io/badge/Biomni-A1%20Agent-4B8BBE)
![DepMap](https://img.shields.io/badge/Data-DepMap%2FCCLE-E4405F)
![SynLethDB](https://img.shields.io/badge/Data-SynLethDB%202.0-6E40C9)
![cBioPortal](https://img.shields.io/badge/API-cBioPortal-1F6FEB)
![NCBI](https://img.shields.io/badge/API-NCBI%20Entrez-336791)
![STRING](https://img.shields.io/badge/API-STRING--DB-05998B)

| 구분 | 사용 기술 |
|---|---|
| 통계 | Welch t-test(단측/양측), Benjamini–Hochberg FDR, Cohen's d, permutation test, 4-parameter logistic 회귀 |
| 데이터 | DepMap CRISPR gene effect(Chronos)·gene dependency 확률·Model·발현(log2 TPM+1), SynLethDB 2.0, Broad Repurposing Hub |
| 외부 API | cBioPortal(CCLE 변이/CN 콜), NCBI Entrez E-utilities, STRING-DB v12 |
| 시각화 | matplotlib/seaborn(PNG), 의존성 없는 인라인 SVG(HTML 보고서) |
| 에이전트 | Biomni A1(LangGraph), `module2api` 도구 레지스트리 |
| 품질 관리 | Ruff, 모든 출력에 provenance·QC 경고 동반 |

---

## 프로젝트 구조

```
biomni/
├── tool/
│   ├── synthetic_lethality.py        # 발굴·검증·반증 + 난소암 층화 + 시각화 (8개 도구)
│   ├── organoid_sl.py                # 오가노이드 관점 발굴·전이성 평가·실험 설계 (6개 도구)
│   ├── sl_multichannel.py            # DAISY식 다채널 랭킹, FDR 게이트 없음 (1개 도구)
│   ├── pdac_translation.py           # 췌장암 변환연구: 드라이버→약물→PDO→채점→보고서 (6개 도구)
│   └── tool_description/             # 각 모듈의 에이전트용 스키마 (JSON 형태)
├── utils.py                          # read_module2api() 필드 목록에 4개 모듈 등록
run_pdac_agent.py                     # PDAC 8단계 파이프라인 (예측 1-6 / 실험 채점 7-8 / 에이전트)
serve_report.py                       # 리포트 서버 (파이프라인과 분리, 표준 라이브러리만 사용)
run_sl_agent.py                       # 단계별 실행 / A1 자연어 실행
run_agent.py                          # 췌장암 KRAS 세포주·오가노이드 관점 병렬 산출
run_ferroptosis_agent.py              # 워크플로를 지정하지 않고 에이전트가 직접 계획하는 예시
data/biomni_data/data_lake/           # DepMap·SynLethDB·Repurposing Hub 스냅샷
pdac_run/                             # 파이프라인 출력 (CSV + HTML 보고서)
```

도구 등록은 Biomni 규약을 그대로 따릅니다. `biomni/utils.py`의 `read_module2api()` 필드 목록에
`organoid_sl`, `pdac_translation`, `sl_multichannel`, `synthetic_lethality`가 올라가 있어
`A1` 초기화 시 **21개 도구**가 자동으로 레지스트리에 등록됩니다.

---

## 도구 목록

### `synthetic_lethality` — 발굴, 검증, 반증 (8개)

| 도구 | 역할 |
|---|---|
| `discover_synthetic_lethal_candidates` | 변이/야생형 층화 후 전 유전자 Welch t-test, pan-essential 필터 |
| `validate_sl_candidates_with_pubmed` | Entrez esearch+efetch, 0–100 문헌 지지 점수 (반증 전용 질의 계층 포함) |
| `analyze_ppi_network_for_sl` | STRING 직접 edge + 공유 파트너 + 기능 enrichment |
| `check_dependency_confounders` | **반증 전용** — 파라로그 손실·발현 바이오마커·lineage 교란 검사 |
| `generate_sl_evidence_dossier` | 위 단계 통합, Go/Hold/No-go 판정 |
| `stratify_ovarian_cancer_dependency_by_mutation` | 난소암 층화 — 변이/증폭/결손 3가지 유전형 이벤트, SynLethDB 주석 |
| `plot_dependency_boxplot` | 그룹별 세포주 의존성 박스+스트립 (패널별 Δ·p 재계산) |
| `plot_dependency_volcano` | 효과 차이 vs −log10(p), pan-essential 별도 마커, SynLethDB 파트너 표시 |

### `organoid_sl` — 오가노이드 관점 (6개)

`discover_subtype_masked_sl_candidates` · `discover_serum_masked_dependencies` ·
`discover_allele_resolved_sl_candidates`(KRAS G12D/G12V 등 allele별) · `assess_organoid_transferability` ·
`design_organoid_sl_experiment` · `rank_organoid_sl_candidates`

### `sl_multichannel` — 다채널 랭킹 (1개)

`discover_sl_multichannel` — essentiality·co-expression·co-inactivation 세 채널을 랭크 결합.
소규모 코호트에서 genome-wide FDR이 진짜 상호작용까지 버리는 문제를 피하기 위해 **임계값이 아닌 순위**로 보고합니다.

### `pdac_translation` — 췌장암 변환연구 (9개)

| 도구 | 역할 |
|---|---|
| `profile_pdac_driver_landscape` | 각 드라이버가 **검정 가능한지** 먼저 판정 (TESTABLE / UNDERPOWERED / UNTESTABLE) |
| `discover_sl_pan_cancer_with_context` | **전략 A** — pan-cancer로 검정력 확보 후 PDAC 맥락에서 존재 여부 재확인. allele 특이(G12D), 우세 계통 제거 재검정 포함 |
| `discover_comutation_stratified_sl` | **전략 B** — KRAS 변이주를 TP53/SMAD4/CDKN2A 상태로 층화. 드라이버가 아니라 **공변이 유전형**의 의존성 |
| `compare_allele_specific_dependencies` | G12D vs G12V — 세 대비(A/WT, B/WT, **A/B 직접**)로 공통 파트너와 allele 특이 파트너 구분, MRTX1133 병용 축 주석 |
| `map_sl_candidates_to_drugs` | 세 소스(큐레이션 RAS 세트 + Repurposing Hub + ChEMBL 실시간) 매핑. MOA가 저해 기전인 화합물만 약물 arm으로 집계 |
| `analyze_organoid_drug_response` | 측정된 PDO 용량반응 4PL 피팅 → IC50/AUC/Emax, 유전형군 비교 |
| `analyze_crispr_validation` | 오가노이드 CRISPR KO 생존율 → log2FC + 유전형 선택성 판정 |
| `compare_prediction_with_experiment` | **예측 채점** — 반증된 예측 우선 보고, precision@k + permutation null |
| `generate_pdac_report` | 위 결과를 차트·표·실제 PubMed 참고문헌이 담긴 HTML 보고서로 생성 |

---

## 설계 원칙

1. **수치 판단은 코드가, 서술은 LLM이** — 통계·임계값·QC는 명시적 규칙으로 계산하고 LLM은 계획과 설명만 담당합니다.
2. **모든 실행에 provenance** — 데이터 파일·스냅샷 날짜·표본 수·API 엔드포인트·통계 가정을 함께 반환합니다.
3. **신규성과 반증을 구분** — 문헌이 없는 후보는 중립에서 시작하고, **명시적 반증 문헌만** 감점합니다.
4. **능동적 반증 탐색** — PubMed 질의를 질환 특이 / SL 특이 / **반증 특이** 3계층으로 나눕니다.
5. **빈 결과도 결과** — FDR을 통과한 유전자가 없으면 빈 표 대신 runner-up과 **"컷오프를 올리지 말고 검정 대상을 줄여라"**는
   안내를 출력합니다. 소규모 코호트에서 BH 보정을 통과하지 못하는 것은 검정력의 한계이지 효과가 없다는 근거가 아닙니다.
6. **마지막 단계는 "아니오"라고 말할 수 있어야 한다** — 채점 도구는 확인된 예측보다 **반증된 예측을 먼저** 출력합니다.

---

## 설치 및 실행

### 1. 설치

```bash
git clone https://github.com/Dandanzzi/Biomni.git
cd Biomni
pip install -e .
```

시각화에 `matplotlib`, `seaborn`이 필요하며 둘 다 `requirements.txt`에 포함돼 있습니다.
parquet 읽기에 `pyarrow`가 필요합니다(시스템 파이썬에 없으면 아래처럼 venv를 쓰세요).

> 이 저장소에서는 `venv/bin/python`을 사용하며, **반드시 저장소 루트에서 실행**합니다
> (데이터 레이크 경로가 `./data/biomni_data/data_lake` 상대경로이기 때문입니다).
> 다른 위치에서 돌리려면 각 도구의 `data_lake_path` 인자에 절대경로를 주세요.

### 2. 데이터 준비

| 파일 | 용도 | 필수 여부 |
|---|---|---|
| `DepMap_CRISPRGeneEffect.csv` | 유전자 의존성 (Chronos) | 필수 |
| `DepMap_Model.csv` | 세포주 계통/암종 주석 | 필수 |
| `DepMap_CRISPRGeneDependency.csv` | 의존 확률(0–1), essentiality floor | `sl_multichannel` |
| `DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv` | 발현 바이오마커 상관 | 반증·전이성 평가 |
| `synlethdb_human_sl.parquet` | SynLethDB 2.0 인간 SL 쌍 37,943개 | SynLethDB 주석 |
| `broad_repurposing_hub_phase_moa_target_info.parquet` | 유전자→약물·임상단계·MOA | `map_sl_candidates_to_drugs` |

> **Repurposing Hub 스냅샷은 KRAS 저해제 시대 이전입니다** — sotorasib·adagrasib·MRTX 계열이 0건이고 SOS1 주석도 없습니다.
> 그래서 `map_sl_candidates_to_drugs`는 세 소스를 합칩니다: 모듈에 포함된 **큐레이션 RAS 경로 약물표**(KRAS G12C/G12D,
> SOS1, SHP2, WRN, TEAD 등, 날짜 명시), Hub, 그리고 **ChEMBL 실시간 조회**. ChEMBL이 응답하지 않으면 오프라인 소스로
> 저하되며 그 사실을 QC 경고와 provenance에 표시합니다 — "약물 없음"과 "조회 실패"를 절대 섞지 않습니다.
| `DepMap_OmicsCNGene.csv` | 연속 copy-number | 선택 (없으면 cBioPortal 사용) |

변이와 copy-number 상태는 별도 파일 없이 **cBioPortal CCLE API**에서 자동 조회합니다.
로컬에 DepMap 변이/CN 파일이 있으면 그쪽을 우선 사용하며, `mutation_csv_path`로 직접 지정할 수도 있습니다
(변이: `ModelID, HugoSymbol[, ProteinChange]` / CN: `ModelID, HugoSymbol[, Alteration]`).

> SynLethDB 스냅샷의 `score` 컬럼은 전 행이 1.0이라 신뢰도 가중에 쓸 수 없습니다.
> 도구는 `source`를 근거 등급(실험 > 대규모 스크린 > text mining > 계산 예측)으로 매핑해 사용하며, 이 사실을 QC 경고에 표시합니다.

인터넷이 필요한 구간: cBioPortal(변이/CN), NCBI Entrez(문헌), STRING(네트워크).

### 3. PDAC 파이프라인 실행

연구 질문: *췌장암의 유전적 취약성과 합성치사 관계를 발굴하고, 예측된 약물 효과를 환자유래 오가노이드에서 기능적으로 검증할 수 있는가?*

| 단계 | 도구 | 산출 |
|---|---|---|
| 1 | `profile_pdac_driver_landscape` | 어떤 드라이버가 검정 가능한가 |
| 2 | `discover_synthetic_lethal_candidates` | DepMap 유전형 대비 |
| 3 | `discover_sl_multichannel` | 순위 기반 채널 (`--multichannel`) |
| 4 | `rank_organoid_sl_candidates` | 오가노이드를 쓸 가치가 있는 후보 |
| 5 | `map_sl_candidates_to_drugs` | 약물 arm vs CRISPR arm 분기 |
| 6 | `design_organoid_sl_experiment` | 반증 기준이 포함된 프로토콜 |
| 7 | `analyze_organoid_drug_response` / `analyze_crispr_validation` | 측정 데이터 분석 |
| 8 | `compare_prediction_with_experiment` | 예측 채점 |

**예측 단계 (LLM 불필요, API 비용 없음, 약 2–3분)**

```bash
venv/bin/python run_pdac_agent.py --mode direct --stages 1-6 --output-dir ./pdac_run
```

**실험 데이터 채점 (측정 후)**

```bash
venv/bin/python run_pdac_agent.py --mode direct --stages 7-8 \
  --crispr-csv crispr_validation.csv --genotype-csv pdo_genotype.csv \
  --screen-csv pdo_screen.csv --prediction-csv ./pdac_run/prediction_depmap.csv --top-k 5
```

입력 CSV 형식(컬럼명은 모두 인자로 바꿀 수 있습니다):

```
# pdo_screen.csv          organoid,drug,concentration_um,viability,replicate
# crispr_validation.csv   organoid,gene,viability,replicate   (비표적 대조는 gene=NTC)
# pdo_genotype.csv        organoid,genotype                   (예: KRAS-G12D / KRAS-WT)
```

**HTML 보고서**

보고서는 위 두 명령이 끝날 때 **자동으로** `pdac_run/pdac_report.html`에 생성됩니다
(`--no-report`로 생략, `--no-literature`로 PubMed 조회 생략, `--email`로 Entrez 연락처 지정).

보고서만 다시 만들려면:

```bash
venv/bin/python -c "
from biomni.tool.pdac_translation import generate_pdac_report
print(generate_pdac_report(run_dir='./pdac_run', email='you@example.org'))
"
```

### 5. 리포트 서버

생성된 HTML을 브라우저로 보려면 `serve_report.py`를 띄웁니다. **파이프라인과 분리된 독립 스크립트**이고
표준 라이브러리만 쓰므로 biomni 패키지도, 외부 서비스 로그인도 필요하지 않습니다.

```bash
python serve_report.py                    # ./pdac_run 을 http://localhost:8000/ 에 서비스
python serve_report.py --port 8080        # 포트 변경
python serve_report.py --host 0.0.0.0     # 같은 네트워크의 다른 PC에서도 접근
python serve_report.py --background       # 백그라운드 실행 후 프롬프트 복귀
python serve_report.py --status           # 실행 여부 확인
python serve_report.py --stop             # 종료
```

`/` 로 접속하면 디렉토리의 리포트·그림·CSV 목록이 나오고, 리포트 파일로 바로 들어갈 수 있습니다.
VS Code 원격(SSH) 세션이면 **포트가 자동 전달**되므로 PORTS 탭의 지구본 아이콘으로 로컬 브라우저에서 열립니다.

> `--host 0.0.0.0`은 이 서버에 접근 가능한 모든 사람에게 리포트를 공개합니다. 기본값은 `127.0.0.1`입니다.

드라이버 검정 가능성·후보 선택성·공략 가능성·문헌 지지도를 인라인 SVG 차트(호버 툴팁 + 원본 데이터 표)로 그리고,
각 후보의 **실제 PubMed 참고문헌을 PMID 링크로** 싣습니다. 외부 스크립트를 불러오지 않아 오프라인에서도 열립니다.
실험 데이터가 아직 없으면 해당 섹션을 "대기"로 표시하고 사전 지정 기준만 명시합니다.

### 6. 그 밖의 실행 경로

```bash
venv/bin/python run_sl_agent.py --mode direct          # 범용 SL 체인, 도구별 원시 출력
venv/bin/python run_agent.py                           # 세포주 관점 + 오가노이드 관점 병렬 산출
venv/bin/python run_pdac_agent.py --mode agent --guided # A1이 자연어 질문에서 직접 도구 체인 계획 (LLM 키 필요)
```

```python
from biomni.agent import A1

agent = A1(path="./data", llm="claude-sonnet-4-5-20250929")
agent.go("췌장암에서 KRAS 합성치사 후보를 찾고 오가노이드 검증 실험을 설계해 줘")
```

---

## 실행 예시

### 1) 발굴보다 먼저: 이 드라이버가 검정 가능한가

`profile_pdac_driver_landscape` 실제 출력(DepMap 스냅샷, 췌장암 세포주 47종):

| driver | alteration | 변이군 | 대조군 | 빈도 | 판정 |
|---|---|---|---|---|---|
| TP53 | mutation | 36 | 8 | 82% | TESTABLE |
| SMAD4 | deletion | 26 | 10 | 72% | TESTABLE |
| CDKN2A | deletion | 23 | 13 | 64% | TESTABLE |
| **KRAS** | **mutation** | **40** | **4** | **91%** | **UNDERPOWERED** |
| BRCA2 | mutation | 2 | 42 | 5% | UNTESTABLE |

가장 유명한 드라이버가 바로 그 빈도 때문에 **야생형 대조군이 4종뿐**입니다. 발굴을 돌리기 전에 이 사실을 알려주는 것이
이 도구의 목적이며, 대안(SMAD4 결손·CDKN2A 결손, 또는 allele별/순위 기반 분석)도 함께 제시합니다.

### 2) 통계 1순위가 반증으로 탈락한다

KRAS 분석 상위 후보 중 **VPS4A(q=0.005)는 `check_dependency_confounders`에서 탈락**합니다.

```
PARALOG ALERT: VPS4B expression (r=+0.315, rank 6 of 19199) predicts this dependency
far better than KRAS (rank 16167) - this looks like a paralog-loss dependency.

VERDICT: CONFOUNDED - the dependency tracks VPS4B expression, not KRAS status.
```

같은 후보에 대해 보고서의 문헌 섹션도 독립적으로 같은 방향을 가리킵니다 — VPS4A는 검색된 논문 2편 모두가
KRAS와 함께 다룬 논문이 아닌 **단순 동시 언급**이라 지지 점수 4/100입니다.

**문헌 검증에서도 같은 원리가 작동합니다.** 반증 전용 질의 계층을 추가하기 전 STK33은 89/100 "WELL SUPPORTED"였지만,
추가 후 [PMID 21742770](https://pubmed.ncbi.nlm.nih.gov/21742770/) (*STK33 kinase activity is nonessential in
KRAS-dependent cancer cells*)을 회수해 **SUPPORTED BUT CONTESTED**로 강등됩니다.

### 3) 후보 9개 중 4개가 오가노이드 투입 대상에서 제외된다

`rank_organoid_sl_candidates`는 각 후보에 클래스를 부여하고, 상위 랭크를 받을 수 없는 클래스를 분리합니다.
췌장암 KRAS 실행(FDR 통과 9개) 결과:

| 유전자 | 클래스 | 이유 |
|---|---|---|
| TEAD1 | SL-CANDIDATE | 변이주 32% 의존 / 야생형 0%, 전체 17% — 선택적 |
| DHFR, ADAR | NO-WINDOW | 야생형·전체 세포주에서도 의존(DHFR 전체 72%, ADAR 51%) |
| VPS4A, ARF4 | PARALOG-CONFOUNDED | 의존성이 driver가 아니라 파라로그 발현을 따라감 |

통계 2순위였던 **VPS4A가 여기서도 PARALOG-CONFOUNDED로 걸립니다** — 위 2)의 반증 검사와 독립적으로 같은 결론에 도달합니다.
랭킹에 남은 5개(TEAD1, SREBF1, NDE1, FADD, CHKA) 중 최종 픽 TEAD1은 Repurposing Hub에 저해제가 없어
**약물 arm이 아니라 CRISPR arm**으로 라우팅됩니다.

### 4) KRAS 이분법을 우회하는 두 가지 층화

PDAC에서 KRAS 야생형은 4종뿐이라 어떤 통계로도 대비가 생기지 않습니다. 두 가지 우회로가 있습니다.

**전략 A — pan-cancer 후 맥락 확인** (`discover_sl_pan_cancer_with_context`). KRAS G12D를 전 암종에서 모으면
**47 vs 791**이 되어 q<0.25 통과 유전자가 0개에서 **537개**로 늘어납니다. 양성 대조도 복원됩니다 — KRAS 자신이
delta −1.516, q=0.0000으로 1위이고 PDAC G12D 세포주 100%가 의존합니다.

| gene | MUT | WT | delta | q | 계통제거 | PDAC 평균 | PDAC 의존 | 판정 |
|---|---|---|---|---|---|---|---|---|
| KRAS | −1.899 | −0.383 | −1.516 | 0.0000 | −1.274 | −2.257 | 100% | PRESENT |
| KLF5 | −0.721 | −0.345 | −0.376 | 0.0449 | −0.365 | −0.744 | 58% | PRESENT |
| RAB10 | −0.781 | −0.413 | −0.367 | 0.0133 | −0.334 | −0.832 | 79% | PRESENT |
| CCND1 | −1.261 | −0.925 | −0.337 | 0.1236 | **−0.193** | −1.476 | 100% | PRESENT |
| TCF7L2 | −0.341 | −0.063 | −0.279 | 0.0449 | −0.333 | −0.261 | 21% | **ABSENT** |

'계통제거' 열은 변이군의 우세 계통(췌장, 40%)을 빼고 같은 대비를 다시 돌린 값입니다. CCND1처럼 이 값이 0 쪽으로
무너지면 유전형 효과가 아니라 조직 효과였다는 뜻입니다. TCF7L2는 pan-cancer에서 유의하지만 PDAC 세포주에서는
의존성 자체가 없어 ABSENT로 걸립니다.

**전략 B — 공변이 층화** (`discover_comutation_stratified_sl`). 모든 세포주가 KRAS 변이인 상태에서 두 번째 변이로 나눕니다.

| 층화 | 군 크기 | q<0.25 | 최종 후보 |
|---|---|---|---|
| TP53 변이 vs 야생형 | 34 vs 6 | 10 | IRS2 |
| SMAD4 결손 vs 정상 | 25 vs 9 | 5 | VPS4A |
| CDKN2A 결손 vs 정상 | 22 vs 12 | 0 | — |
| SMAD4 변이 vs 야생형 | 11 vs 29 | 0 | — |

SMAD4 결손군의 유일한 후보 VPS4A는 **교란입니다.** SMAD4는 18q21.2, 그 파라로그 VPS4B는 18q21.33에 있어
18q 대규모 결손이면 함께 사라집니다. 실제로 확인됩니다 — SMAD4 결손 세포주의 VPS4B 발현이 유의하게 낮습니다
(4.41 vs 5.11 log2 TPM+1, Welch p=0.0008). 즉 "SMAD4 결손 합성치사"가 아니라 **동반 결손된 VPS4B의 파라로그 의존성**입니다.

### 5) KRAS allele은 하나의 유전형이 아니다

`compare_allele_specific_dependencies`는 **세 가지 대비**를 돌립니다. A vs 야생형, B vs 야생형만 보면 한쪽이
단지 표본이 커서 통과했을 수 있으므로, **A vs B 직접 대비**가 확인해 줄 때만 allele 특이로 분류합니다.

pan-cancer KRAS G12D 47 / G12V 31 / 야생형 791 기준 결과:

| 분류 | 유전자 | 근거 |
|---|---|---|
| **G12D 특이** | KLF5, DDX39B | 야생형과 G12V 둘 다보다 강함 (직접 대비 p=0.042 / 0.014) |
| G12D 단독(미확인) | CCND1, CTNNB1, DOCK5, WNK1 등 10개 | 야생형 대비는 통과하나 G12V와 분리되지 않음 |
| 공통 | **KRAS**, RAB10, CFLAR, UROD, TUBB4B | allele 구분 없이 KRAS 변이 전반의 파트너 |
| G12V 특이 | TRPM7, SCAP, RAB6A | |

양성 대조가 정확합니다 — KRAS 자신이 **공통**으로 분류되고(두 allele 모두 q=0.0000) 직접 대비는 유의하지
않습니다(dA−B=−0.126, p=0.24). 두 allele 모두 KRAS에 똑같이 의존하므로 이게 맞는 답입니다.

### 6) 난소암 층화 — 변이로는 보이지 않는 드라이버

`stratify_ovarian_cancer_dependency_by_mutation`은 변이뿐 아니라 **증폭·결손**으로도 층화합니다.
HGSOC의 핵심 드라이버 CCNE1은 증폭 이벤트라 변이 콜로는 군이 만들어지지 않습니다.

```
CCNE1-AMPLIFIED : NIHOVCAR3, ONCODG1, COV318, KURAMOCHI, JHOS4  (n=5)
CCNE1-NEUTRAL   : CAOV4, OAW28, JHOS2, COV362, TYKNU, HEYA8, ...  (n=9)
Self-dependency control (WEAK): CCNE1 delta=-0.764, p1=5.0e-02
```

---

## 알려진 한계

이 도구는 가설 생성기이며, 아래 한계를 리포트에 명시적으로 출력합니다.

- **단일 perturbation 상관 근거**입니다. 이중 녹아웃·isogenic 검증이 아니므로 확정된 SL 상호작용이 아닙니다.
- **소규모 코호트에서는 genome-wide FDR을 통과하는 것이 원리적으로 어렵습니다.** 난소암 59종 코호트에서 BH 게이트를
  통과하려면 Δ < −0.44가 필요한데 실제 타겟은 −0.2 ~ −0.4 수준입니다. 유전자셋을 좁히거나 순위 기반 채널을 쓰세요.
- 단일 DepMap release에서 파생된 결과는 **독립 증거로 중복 계산하지 않습니다**(Sanger Project Score 등 외부 교차검증 권장).
- 문헌 점수는 **문서화 정도이지 진위가 아닙니다.** 보고서는 드라이버와 후보를 실제로 함께 다룬 논문 수를 따로 표시하므로,
  점수가 높아도 `pair=0`이면 단순 동시 언급입니다.
- 변이 콜은 기능 결손과 VUS를, 단일대립과 양대립 손실을 구분하지 않습니다. BRCA1 프로모터 메틸화처럼 변이가 아닌
  기전은 야생형군에 섞여 대비를 희석시킵니다.
- 세포주·오가노이드 의존성은 환자에서의 therapeutic window를 보장하지 않습니다.
- **채점 단계를 반증 가능하게 만들려면 실험 설계가 받쳐줘야 합니다.** 예측 상위만 검증하면 precision@k가 base rate와
  같아져 우연과 구분되지 않습니다. KO 패널에 **예측이 낮게 평가한 유전자도 포함**하세요.

---

## 참고

- Huang K, et al. *Biomni: A General-Purpose Biomedical AI Agent.* bioRxiv (2025)
- DepMap Consortium — [Cancer Dependency Map portal](https://depmap.org/portal/)
- SynLethDB 2.0 — [인간 합성치사 쌍 데이터베이스](https://synlethdb.sist.shanghaitech.edu.cn/)
- Broad Institute — [Drug Repurposing Hub](https://clue.io/repurposing)
- Jerby-Arnon L, et al. *Predicting cancer-specific vulnerability via data-driven detection of synthetic lethality.* Cell (2014)
- Neggers JE, et al. *Synthetic lethal interaction between the ESCRT paralog enzymes VPS4A and VPS4B.* Cell Reports (2020)
