# Biomni Synthetic Lethality Toolkit

**암 맥락별 합성치사(Synthetic Lethality) 후보 발굴부터 반증 검사까지 — Biomni 호환 도구 생태계**

DepMap CRISPR 의존성 데이터로 후보를 발굴하고, PubMed 문헌과 STRING 단백질 네트워크로 교차 검증한 뒤,
**"이 후보가 정말 driver 변이 때문인가"를 능동적으로 반증**하여 실험 가능한 evidence dossier를 생성합니다.

> 이 도구는 "점수를 잘 내는 예측기"가 아니라 **"어떤 후보를 왜 먼저 검증해야 하는가"에 답하는 의사결정 계층**을 목표로 합니다.
> 그래서 지지 근거만 모으지 않고, pan-essentiality·파라로그 손실·lineage 교란·문헌상 반증을 적극적으로 찾아 후보를 탈락시킵니다.

---

## 기술 스택

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
![pandas](https://img.shields.io/badge/pandas-2.x-150458?logo=pandas&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-2.x-013243?logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-1.15-8CAAE6?logo=scipy&logoColor=white)
![Requests](https://img.shields.io/badge/Requests-HTTP-2C5BB4)
![Ruff](https://img.shields.io/badge/lint-Ruff-D7FF64?logo=ruff&logoColor=black)

![Biomni](https://img.shields.io/badge/Biomni-A1%20Agent-4B8BBE)
![DepMap](https://img.shields.io/badge/Data-DepMap%2FCCLE-E4405F)
![cBioPortal](https://img.shields.io/badge/API-cBioPortal-1F6FEB)
![NCBI](https://img.shields.io/badge/API-NCBI%20Entrez-336791)
![STRING](https://img.shields.io/badge/API-STRING--DB-05998B)

| 구분 | 사용 기술 |
|---|---|
| 통계 분석 | Welch t-test, Benjamini–Hochberg FDR, Cohen's d, 벡터화 Pearson 상관, 아형 z-score |
| 로컬 데이터 | DepMap CRISPR gene effect / Model / 발현, GTEx 정상조직, BioGRID 유전상호작용, Ensembl gene_info |
| 외부 API | cBioPortal (CCLE 변이 콜), NCBI Entrez E-utilities (esearch/efetch), STRING-DB v12 |
| 에이전트 | Biomni A1 (LangGraph 기반), `module2api` 도구 레지스트리 (총 232개 중 8개가 이 툴킷) |
| 품질 관리 | Ruff (lint + format), 명시적 provenance 기록, 결정론적 파이프라인 |

---

## 사용 데이터셋

모든 수치 근거는 아래 6개 로컬 파일과 3개 외부 API에서만 나옵니다. LLM의 사전지식은 수치 판단에 쓰이지 않습니다.

### 로컬 데이터 (`data/biomni_data/data_lake/`)

| 데이터셋 | 규모 | 쓰는 도구 | 무엇을 위해 |
|---|---|---|---|
| `DepMap_CRISPRGeneEffect.csv` | 1,183 모델 × 17,916 유전자 (Chronos 보정) | 전 도구 | 유전자 녹아웃 시 생존 영향. 0=무영향, −1=일반 필수유전자 수준 |
| `DepMap_Model.csv` | 2,116 모델 | 전 도구 | 암종·계통(OncotreeLineage)·배양형태(GrowthPattern)·ModelType 주석 |
| `DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv` | 1,684 모델 × 19,205 유전자 | 교란검사, 아형 층화 | 발현 바이오마커 상관, classical/basal 아형 점수 |
| `gtex_tissue_gene_tpm.parquet` | 54개 정상조직 | `assess_organoid_transferability` | 정상 췌장 발현 = 치료 window 사전 판정 |
| `genetic_interaction.parquet` | 사람 5,557쌍 (1,336개 유전자) | `assess_organoid_transferability` | BioGRID 직교 유전상호작용 근거 |
| `gene_info.parquet` | Ensembl↔심볼 | 내부 | BioGRID의 Ensembl ID를 HUGO 심볼로 매핑 |

> **BioGRID 커버리지 주의**: 사람 유전상호작용은 소수의 조합 CRISPR 스크린에서 온 1,336개 유전자만 덮습니다.
> 대부분의 후보는 "상호작용 없음"이 아니라 **"검사된 적 없음"**이며, 리포트가 그렇게 구분해 표기합니다.

### 외부 API (실행 시 조회)

| API | 용도 | 쓰는 도구 |
|---|---|---|
| cBioPortal `ccle_broad_2019` | 세포주 변이 콜 (KRAS 변이 235건 / 1,739개 프로파일 세포주) | `discover_synthetic_lethal_candidates` |
| NCBI Entrez E-utilities | esearch + efetch XML 초록 파싱 | `validate_sl_candidates_with_pubmed` |
| STRING-DB v12 | 단백질 상호작용, 공유 파트너, 기능 enrichment | `analyze_ppi_network_for_sl` |

변이 상태는 "변이 있음 / 야생형 / **미프로파일**"을 구분합니다. 프로파일되지 않은 세포주를 야생형으로
간주하지 않고 제외하며, 이 때문에 췌장암 47종 중 3종이 분석에서 빠집니다.

---

## 프로젝트 구조

```
biomni/
├── tool/
│   ├── synthetic_lethality.py                  # 세포주 층 (5개 도구 + 내부 헬퍼)
│   ├── organoid_sl.py                          # 오가노이드 층 (3개 도구)
│   └── tool_description/
│       ├── synthetic_lethality.py              # 에이전트용 도구 스키마 (JSON 형태)
│       └── organoid_sl.py
├── utils.py                                    # read_module2api() 필드 목록에 모듈 등록
docs/
└── synthetic_lethality/
    └── README.md                               # 이 문서
run_agent.py                                    # 췌장암 KRAS 원스톱 실행 (후보 → 문헌 → 반증 → 실험계획)
run_sl_agent.py                                 # 단계별 실행 / A1 에이전트 자연어 실행
data/biomni_data/data_lake/                     # DepMap 스냅샷 (아래 "데이터 준비" 참조)
```

도구 등록은 Biomni 규약을 그대로 따릅니다. `biomni/utils.py`의 `read_module2api()` 필드 목록에
`"synthetic_lethality"`와 `"organoid_sl"` 두 줄이 추가되어 있어, `A1` 초기화 시 8개 도구가
자동으로 레지스트리(총 232개)에 올라갑니다. 실제로 Biomni의 tool retriever는 KRAS를 언급하지 않은
자유형 한국어 질문("췌장암과 합성치사 관계에 있는 것을 알려줘")에서 이 도구들을 선택했습니다.

---

## 프로그램 흐름

`run_agent.py`는 LLM을 쓰지 않는 **결정론적 파이프라인**입니다. 같은 입력이면 항상 같은 결과가 나오고,
"왜 이 후보를 골랐는가"가 전부 코드에 고정된 규칙으로 설명됩니다.

```
사용자 질문: "췌장암에서 KRAS 변이 기반 합성치사 후보를 찾아줘"
      │
[1] 후보 발굴  discover_synthetic_lethal_candidates
      │  DepMap 췌장암 47종 → cBioPortal로 KRAS 변이 층화 (변이 40 / 야생형 4 / 미프로파일 3)
      │  17,787개 유전자 Welch t-test → BH-FDR → pan-essential 제거
      │  양성대조: KRAS 자기 자신이 선택적 의존성으로 회수되는지 확인
      ↓  9개 후보
[2] 해석 + 우선순위 선별   (규칙: pan-essential<40%, q<0.1, 변이의존≥30%)
      │  각 유전자가 왜 통과/탈락했는지 근거와 함께 화면 출력
      ↓  5개
[3] 문헌 검증  validate_sl_candidates_with_pubmed
      │  Entrez 3계층 질의: 질환특이 / SL특이 / **반증특이**
      ↓  0–100 지지 점수 + 반증 문헌 PMID
[4] 반증 검사  check_dependency_confounders          ← 실험비를 쓰기 전의 관문
      │  driver와의 co-dependency / 발현 바이오마커 / 파라로그 스캔 / lineage 집중도
      ↓  CONFOUNDED 후보는 여기서 탈락 (예: VPS4A → VPS4B 파라로그)
[5] 오가노이드 관점
      │  (5a) discover_subtype_masked_sl_candidates — 아형 평균화로 상쇄된 후보 발굴
      │  (5b) assess_organoid_transferability      — 오가노이드 배지에서 가려질지 판정
      ↓  통합 요약표: 세포주 관점 vs 오가노이드 관점
[6] 실험 설계
         print_experiment_plan            — 2D isogenic 검증 프로토콜
         design_organoid_sl_experiment    — PDO 프로토콜 (배지 수정·정상 대조군·반증 기준)
```

각 단계는 파일이 아니라 **research log 문자열**을 화면에 출력합니다. 이것이 Biomni의 도구 규약이며,
에이전트가 그대로 읽고 다음 단계를 판단할 수 있는 형태입니다.

---

## 핵심 기능

### 8개 도구

**세포주 층** (`biomni/tool/synthetic_lethality.py`)

| 도구 | 역할 | 주요 출력 |
|---|---|---|
| `discover_synthetic_lethal_candidates` | DepMap 세포주를 변이/야생형으로 층화해 전 유전자 Welch t-test | 효과크기, p/q값, 선택도, pan-essential 지표, QC 경고 |
| `validate_sl_candidates_with_pubmed` | NCBI Entrez esearch + efetch로 초록 XML 파싱 | 논문 목록(PMID/연도/저널/초록), 0–100 문헌 지지 점수 |
| `analyze_ppi_network_for_sl` | STRING-DB 직접 edge + 공유 파트너 + 기능 enrichment | 근거 채널별 점수, 기능적 근접성 분류 |
| `check_dependency_confounders` | **반증 전용** — 교란요인이 의존성을 설명하는지 검사 | co-dependency, 발현 바이오마커, 파라로그 스캔, lineage 집중도 |
| `generate_sl_evidence_dossier` | 위 단계를 통합 | 신뢰등급, 지지/반대 근거, 최소 검증 실험, Go/Hold/No-go |

**오가노이드 층** (`biomni/tool/organoid_sl.py`) — 세포주에서 구조적으로 안 보이는 것을 다룹니다

| 도구 | 역할 | 주요 출력 |
|---|---|---|
| `discover_subtype_masked_sl_candidates` | 아형(classical/basal) 별로 나눠 재검정 — 통합 분석이 상쇄시킨 후보 발굴 | 아형 내에서만 유의한 유전자, 아형별 검정력 |
| `assess_organoid_transferability` | 오가노이드에서도 보일지 예측 (배지 niche 인자, 부착 의존성, 정상조직 window, 직교 유전상호작용) | ORGANOID-ENHANCED / MASKED / WEAKENED / TRANSFERABLE + **배지 수정 권고** |
| `design_organoid_sl_experiment` | PDO 검증 프로토콜 생성 | 모델 패널·배지·CRISPR 전달·3D 판독·검정력·사전 반증 기준 |

### 설계 원칙

0. **결론을 먼저, 10,000자 안에** — Biomni는 도구 출력을 첫 10,000자에서 자릅니다
   (`biomni/agent/a1.py`). 초기 버전은 22,357자를 반환하고 Go/Hold/No-go를 맨 끝에 두어
   **에이전트가 판정을 아예 보지 못했고**, 그 결과 pan-essential 유전자를 "최우선 타겟"으로 서술했습니다.
   현재 `generate_sl_evidence_dossier`는 순위표와 권고를 맨 앞에 두고 9,700자 이내로 반환하며,
   원시 단계 로그는 `include_stage_logs=True`일 때만 덧붙입니다.
1. **수치 판단은 코드가, 서술은 LLM이** — 통계·임계값·QC는 모두 명시적 규칙으로 계산하며 LLM은 계획과 설명만 담당합니다.
2. **모든 실행에 provenance 기록** — 데이터 파일·스냅샷 날짜·행/열 수·표본 수·API 엔드포인트·통계 가정을 함께 반환합니다.
3. **신규성과 반증을 구분** — 문헌이 없는 후보는 중립 기준점에서 시작하고, **명시적 반증 문헌만** 감점합니다.
   ("문헌 없음"은 "반증됨"이 아닙니다.)
4. **능동적 반증 탐색** — PubMed 질의를 3계층(질환 특이 / SL 특이 / **반증 특이**)으로 나눠, 지지 근거뿐 아니라
   비재현·비필수성 보고를 일부러 찾아냅니다.
5. **반증을 LLM의 재량에 맡기지 않음** — 교란요인 검사는 `generate_sl_evidence_dossier` 안에서 자동 실행되며,
   CONFOUNDED 판정은 통계가 아무리 강해도 **강제로 No-go**가 됩니다. 에이전트가 이 단계를 "건너뛰기로 결정"할 수 없습니다.

---

## 설치 및 실행

### 1. 설치

```bash
git clone https://github.com/Dandanzzi/Biomni.git
cd Biomni
pip install -e .
```

추가 의존성은 없습니다 — `pandas`, `numpy`, `scipy`, `requests`는 Biomni 기본 의존성에 포함되어 있습니다.

> **`-e` 플래그가 중요합니다.** 비-editable로 설치된 `biomni` 패키지가 site-packages에 있으면,
> 저장소 루트 밖에서 스크립트를 실행할 때 그쪽이 먼저 로드되어 이 툴킷 모듈을 찾지 못합니다.
> `python -c "import biomni; print(biomni.__file__)"`가 저장소 경로를 가리키는지 확인하세요.

> Biomni 본체는 Python 3.11 이상을 요구하지만, 이 툴킷 모듈 자체는 3.10에서도 동작하도록 작성되어 있습니다.

### 2. 데이터 준비

`data/biomni_data/data_lake/`에 DepMap 스냅샷이 필요합니다. Biomni `A1` 초기화 시 자동 다운로드되며,
[DepMap 포털](https://depmap.org/portal/download/)에서 직접 받아도 됩니다.

| 파일 | 용도 | 필수 여부 |
|---|---|---|
| `DepMap_CRISPRGeneEffect.csv` | 유전자 의존성 (Chronos) | 필수 |
| `DepMap_Model.csv` | 세포주 계통/암종 주석 | 필수 |
| `DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv` | 발현 바이오마커·아형 점수 | 반증 검사와 오가노이드 층에 필요 |
| `gtex_tissue_gene_tpm.parquet` | 정상조직 발현 (치료 window) | 오가노이드 전이성 판정에만 |
| `genetic_interaction.parquet` + `gene_info.parquet` | 직교 유전상호작용 근거 | 선택 (없으면 해당 검사만 생략) |

변이 상태는 별도 파일 없이 **cBioPortal CCLE API**에서 자동 조회합니다.
`DepMap_OmicsSomaticMutations.csv`가 데이터 레이크에 있으면 그쪽을 우선 사용하며,
`mutation_csv_path`로 직접 지정할 수도 있습니다 (`ModelID, HugoSymbol[, ProteinChange]` 컬럼).

인터넷 접속이 필요한 구간: cBioPortal(변이), NCBI Entrez(문헌), STRING(네트워크).
통계 단계만 오프라인으로 돌리려면 `--skip-pubmed`와 로컬 변이 테이블을 사용하세요.

### 3. 실행

**원스톱 실행 (췌장암 KRAS 기본값)**

```bash
python run_agent.py
```

위 [1]~[6] 단계를 한 번에 실행하고 **모든 결과를 화면에 출력합니다.** 파일은 만들지 않습니다.

```bash
python run_agent.py --mutation TP53 --cancer-type "Lung"   # 다른 암종/driver
python run_agent.py --skip-pubmed                          # 문헌 검증 생략 (오프라인)
python run_agent.py --skip-organoid                        # 세포주 관점만
python run_agent.py --csv my_candidates.csv                # 후보 표를 CSV로도 저장 (선택)
```

**단계별 실행 / 에이전트 실행**

```bash
# 각 도구의 원시 출력을 단계별로 확인 (LLM 불필요)
python run_sl_agent.py --mode direct

# 자연어 질문을 A1 에이전트에 그대로 전달 — 도구 선택과 실행 순서는 에이전트가 스스로 계획
python run_sl_agent.py --mode agent --query "췌장암과 합성치사 관계에 있는 것을 알려줘."

# 질문을 생략하면 대화형으로 입력받는다 (빈 줄 입력 시 종료)
python run_sl_agent.py --mode agent
```

에이전트 모드는 워크플로를 지정하지 않습니다. Biomni의 tool retriever가 질문에 맞는 도구를 선택하고,
에이전트가 `<execute>` 파이썬 블록 안에서 도구를 직접 호출합니다.
등록된 도구 목록을 힌트로 덧붙이고 싶으면 `--guided`를 사용하세요.

> **모델 지정**: `--llm` 기본값은 `claude-sonnet-4-5-20250929`입니다.
> Biomni는 LLM 호출 시 `temperature`를 전달하므로, 이 파라미터를 받지 않는 최신 모델에서는 400 오류가 납니다.

**파이썬에서 직접 호출**

```python
from biomni.tool.synthetic_lethality import (
    discover_synthetic_lethal_candidates,
    check_dependency_confounders,
    generate_sl_evidence_dossier,
)

# 1) 후보 발굴
print(discover_synthetic_lethal_candidates("Pancreatic Cancer", "KRAS", top_n=15))

# 2) 실험 전 반증 검사
print(check_dependency_confounders("KRAS", ["VPS4A", "NDE1", "SREBF1"]))

# 3) 전체 파이프라인을 한 번에
print(generate_sl_evidence_dossier("Pancreatic Cancer", "KRAS", top_n=5))
```

**Biomni 에이전트에 장착**

```python
from biomni.agent import A1

agent = A1(path="./data", llm="claude-sonnet-4-5-20250929")
agent.go("췌장암에서 KRAS 변이 기반 합성치사 후보를 찾고 논문으로 검증해 줘")
```

### 어느 실행 방식을 써야 하나

| | `run_agent.py` | `run_sl_agent.py --mode direct` | `run_sl_agent.py --mode agent` |
|---|---|---|---|
| 재현성 | 항상 동일 | 항상 동일 | 매 실행 다름 |
| LLM 개입 | 없음 | 없음 | 계획 + 서술 |
| 후보 선별 | 명시적 규칙 | dossier 통합 점수 | LLM 판단 |
| 용도 | **연구 결과·논문** | 도구별 원시 출력 확인 | Biomni 특성 시연·탐색 |

수치는 세 방식이 동일합니다(같은 도구를 호출하므로). 차이는 신뢰성이며,
에이전트 모드는 실측에서 HOLD 등급 후보를 "최우선 타겟"으로 과장한 적이 있습니다.

---

## 결과 해석

### 출력 컬럼이 뜻하는 것

| 컬럼 | 의미 | 판단 기준 |
|---|---|---|
| `mutant_mean_effect` | 변이 세포주의 평균 gene effect. 0=무영향, −1=일반 필수유전자 수준 | −0.5 아래여야 실제 의존성 |
| `effect_difference` | 변이군 − 야생형군. **음수일수록 합성치사 방향** | −0.3 이하 |
| `cohens_d` | 표준화 효과크기 | \|d\|>0.8이면 큼 |
| `q_value` | 17,787개 유전자 다중검정 보정값 | **q<0.1** (p값이 아니라 이 값을 보세요) |
| `pct_wildtype_dependent` | 야생형 세포주 중 의존 비율 | **0%여야 genotype 선택적** |
| `pct_all_lines_dependent` | 전체 1,183주 중 의존 비율 | 높으면 그냥 필수유전자 → 치료 window 없음 |
| **선택도** | 변이의존% − 전체의존% | 실질적인 우선순위 지표 |

핵심은 마지막 두 컬럼입니다. `effect_difference`가 커도 `pct_all_lines_dependent`가 높으면
정상세포도 함께 죽으므로 치료 표적이 될 수 없습니다.

### 췌장암 KRAS 실행 결과 (DepMap 로컬 스냅샷, 47종 = 변이 40 / 야생형 4 / 미프로파일 3)

| 유전자 | diff | q | %WT 의존 | %전체 의존 | 선택도 | 판정 |
|---|---|---|---|---|---|---|
| **NDE1** | −0.341 | 0.047 | 0% | 30% | **+32.7** | ✅ 최우선 |
| VPS4A | −0.482 | 0.005 | 0% | 20% | +32.9 | ❌ VPS4B 파라로그 교란 |
| **SREBF1** | −0.332 | 0.029 | 0% | 21% | **+21.2** | ✅ 2순위 |
| **TEAD1** | −0.265 | 0.012 | 0% | 17% | **+15.9** | ✅ 오가노이드 최적 |
| DHFR | −0.429 | 0.029 | **50%** | **72%** | +9.3 | ❌ 야생형도 절반 의존 |
| ADAR | −0.496 | 0.122 | 0% | 51% | +8.6 | ❌ FDR 탈락 |
| ARF4 | −0.398 | **0.243** | 25% | 63% | +4.6 | ❌ FDR 탈락 |
| FADD / CHKA | −0.25 / −0.24 | 0.079 / 0.003 | 0% | 12% / 18% | +8.4 / +1.7 | ⚠️ 변이주 20%만 의존 |

**최적 후보군: NDE1, SREBF1, TEAD1.** 특히 TEAD1은 부착 의존성이 압도적(p=9.6e-39)이라
3D 오가노이드에서 효과가 더 크게 나올 유일한 후보입니다.

## 왜 반증 검사가 필요한가

**선택도 1위였던 VPS4A는 `check_dependency_confounders`에서 탈락합니다.**

```
PARALOG ALERT: VPS4B expression (r=+0.315, rank 6 of 19199) predicts this dependency
far better than KRAS (rank 16167) - this looks like a paralog-loss dependency.

VERDICT: CONFOUNDED - the dependency tracks VPS4B expression, not KRAS status.
```

VPS4A 의존성을 설명하는 것은 KRAS가 아니라 파라로그 **VPS4B의 발현**이며(KRAS 발현은 19,199개 중 16,167위),
췌장암 KRAS 야생형 4개 중 3개가 우연히 VPS4B 고발현 세포주였기 때문에 KRAS 효과처럼 보인 것입니다.
반증 검사 없이 진행했다면 KRAS가 아닌 VPS4B 손실 표현형을 실험하게 됩니다.

**문헌 검증에서도 같은 원리가 작동합니다.** 반증 전용 질의 계층을 추가하기 전 STK33은 89/100 "WELL SUPPORTED"였지만,
추가 후 [PMID 21742770](https://pubmed.ncbi.nlm.nih.gov/21742770/) (*STK33 kinase activity is nonessential in
KRAS-dependent cancer cells*)을 회수해 **SUPPORTED BUT CONTESTED**로 강등됩니다.

---

## 오가노이드 층: 세포주가 구조적으로 못 보는 것

세포주 CRISPR 스크린은 유일하게 존재하는 대규모 의존성 자원이지만, 특정 부류의 합성치사를 **원리적으로** 놓칩니다.
이 층은 그 맹점을 계산 가능하게 만들어, 후보를 "그냥 오가노이드도 해보자"가 아니라 **명시된 이유로** PDO 검증에 넘깁니다.

**중요한 한계를 먼저**: 로컬 DepMap에는 췌장 오가노이드 5종(PANFR0069/0233/0368/0402/0420)이 등재되어 있지만
**CRISPR·발현 데이터가 전혀 없습니다.** 이 도구들은 오가노이드 의존성을 측정하지 않습니다.
세포주 근거가 신뢰할 수 없는 지점을 특정하고, 그것을 판정할 오가노이드 실험을 설계합니다.

### 맹점 1: 아형 평균화

췌장암 세포주 47종을 통합 검정하면 classical/basal 아형이 평균화되어 한쪽에만 존재하는 상호작용이 상쇄됩니다.
실제로 돌려보면 더 심각한 사실이 드러납니다:

```
basal-like      mutant=24  wild-type=1     ← KRAS-WT가 1종뿐이라 검정 자체가 불가능
classical-like  mutant=14  wild-type=3
```

**basal-like PDAC에서는 KRAS 합성치사 질문을 기존 세포주로 물을 수 없습니다.** 오가노이드가 필요한 구조적 이유입니다.
classical 아형 내에서만 유의한 후보로 CFLAR, ACLY 등이 나오며, 통합 분석에서는 모두 비유의(p=0.05~0.56)입니다.

### 맹점 2: 배지가 의존성을 가린다

오가노이드 배지는 WNT3A/RSPO1, EGF, FGF10, Noggin, A83-01, Y-27632를 공급합니다.
공급되는 인자의 **상류** 유전자는 녹아웃해도 배지가 구제하므로 **위음성**이 납니다.

```
SMAD4  → ORGANOID-MASKED: A83-01(TGF-β 억제제)이 이미 경로를 차단 → 배지에서 A83-01 제거 필요
PORCN  → ORGANOID-MASKED: 외인성 WNT3A/RSPO1이 리간드 생산 결손을 구제
TEAD1  → ORGANOID-ENHANCED: 부착 배양에서 의존성이 더 강함 (p=9.6e-39) → 3D에서 효과가 더 클 것
```

부착 민감도는 지식이 아니라 **실측**입니다 — DepMap 923개 부착주 vs 162개 부유주 대비로 계산하며,
PTK2(FAK) p=2e-72, ITGB1, YAP1이 예상대로 최상위에 나와 방법이 검증됩니다.

### 맹점 3: 정상 대조군의 부재

세포주 패널에는 짝지어진 정상 조직이 없습니다. 오가노이드에는 있습니다(인접 정상 조직 유래).
GTEx 정상 췌장 발현으로 사전 선별하고(ACLY 25 TPM → WINDOW RISK 플래그),
프로토콜에 **matched normal organoid arm**을 필수 항목으로 넣습니다.

```bash
python -c "
from biomni.tool.organoid_sl import (
    discover_subtype_masked_sl_candidates, assess_organoid_transferability, design_organoid_sl_experiment)
print(discover_subtype_masked_sl_candidates('Pancreatic Cancer', 'KRAS'))
print(assess_organoid_transferability(['TEAD1','ACLY','SMAD4'], 'KRAS'))
print(design_organoid_sl_experiment('TEAD1', 'KRAS'))
"
```

## 알려진 한계

이 도구는 가설 생성기이며, 아래 한계를 리포트에 명시적으로 출력합니다.

- **단일 perturbation 상관 근거**입니다. 이중 녹아웃·isogenic 검증 데이터가 아니므로 확정된 SL 상호작용이 아닙니다.
- **췌장암 KRAS 야생형 세포주는 DepMap 전체에서 4종뿐**이라 p값이 불안정합니다. `cancer_type="pan-cancer"`로 검정력을 높여 재현되는지 확인하세요.
- 단일 DepMap release에서 파생된 결과는 **독립 증거로 중복 계산하지 않습니다** (Sanger Project Score 등 외부 스크린 교차검증 권장).
- 문헌 점수는 **문서화 정도이지 진위가 아닙니다.** 키워드 기반 반증 탐지는 무관한 문맥에서 오탐할 수 있으므로 CONTESTED 판정 시 표시된 PMID를 직접 확인하세요.
- 세포주 의존성은 환자에서의 therapeutic window를 보장하지 않습니다.

---

## 참고

- Huang K, et al. *Biomni: A General-Purpose Biomedical AI Agent.* bioRxiv (2025)
- DepMap Consortium — [Cancer Dependency Map portal](https://depmap.org/portal/)
- Neggers JE, et al. *Synthetic lethal interaction between the ESCRT paralog enzymes VPS4A and VPS4B.* Cell Reports (2020)
- Scholl C, et al. *Synthetic lethal interaction between oncogenic KRAS dependency and STK33 suppression.* Cell (2009) — 및 Babij C, et al. Cancer Research (2011)의 반증 보고
