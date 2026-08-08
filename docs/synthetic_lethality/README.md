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
| 통계 분석 | Welch t-test, Benjamini–Hochberg FDR, Cohen's d, 벡터화 Pearson 상관 |
| 데이터 | DepMap CRISPR gene effect (Chronos), DepMap Model, DepMap 발현 (log2 TPM+1) |
| 외부 API | cBioPortal (CCLE 변이 콜), NCBI Entrez E-utilities (esearch/efetch), STRING-DB v12 |
| 에이전트 | Biomni A1 (LangGraph 기반), `module2api` 도구 레지스트리 |
| 품질 관리 | Ruff (lint + format), 명시적 provenance 기록 |

---

## 프로젝트 구조

```
biomni/
├── tool/
│   ├── synthetic_lethality.py                  # 핵심 구현 (5개 도구 + 내부 헬퍼)
│   └── tool_description/
│       └── synthetic_lethality.py              # 에이전트용 도구 스키마 (JSON 형태)
├── utils.py                                    # read_module2api() 필드 목록에 모듈 등록
docs/
└── synthetic_lethality/
    └── README.md                               # 이 문서
run_agent.py                                    # 췌장암 KRAS 원스톱 실행 (후보 → 문헌 → 반증 → 실험계획)
run_sl_agent.py                                 # 단계별 실행 / A1 에이전트 자연어 실행
data/biomni_data/data_lake/                     # DepMap 스냅샷 (아래 "데이터 준비" 참조)
```

도구 등록은 Biomni 규약을 그대로 따릅니다. `biomni/utils.py`의 `read_module2api()` 필드 목록에
`"synthetic_lethality"` 한 줄이 추가되어 있어, `A1` 에이전트 초기화 시 5개 도구가 자동으로 레지스트리에 올라갑니다.

---

## 핵심 기능

### 5개 도구

| 도구 | 역할 | 주요 출력 |
|---|---|---|
| `discover_synthetic_lethal_candidates` | DepMap 세포주를 변이/야생형으로 층화해 전 유전자 Welch t-test | 효과크기, p/q값, 선택도, pan-essential 지표, QC 경고 |
| `validate_sl_candidates_with_pubmed` | NCBI Entrez esearch + efetch로 초록 XML 파싱 | 논문 목록(PMID/연도/저널/초록), 0–100 문헌 지지 점수 |
| `analyze_ppi_network_for_sl` | STRING-DB 직접 edge + 공유 파트너 + 기능 enrichment | 근거 채널별 점수, 기능적 근접성 분류 |
| `check_dependency_confounders` | **반증 전용** — 교란요인이 의존성을 설명하는지 검사 | co-dependency, 발현 바이오마커, 파라로그 스캔, lineage 집중도 |
| `generate_sl_evidence_dossier` | 위 단계를 통합 | 신뢰등급, 지지/반대 근거, 최소 검증 실험, Go/Hold/No-go |

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
| `DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv` | 발현 바이오마커 상관 | `check_dependency_confounders`에만 필요 |

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

후보 발굴 → 실험 우선순위 선별 → PubMed 검증 → **교란요인 반증 검사** → 최우선 후보의 실험 프로토콜까지 한 번에 출력합니다.

```bash
python run_agent.py --mutation TP53 --cancer-type "Lung"   # 다른 암종/driver
python run_agent.py --skip-pubmed                          # 문헌 검증 생략
python run_agent.py --csv my_candidates.csv                # 후보 표 저장 경로
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

agent = A1(path="./data", llm="claude-sonnet-4-20250514")
agent.go("췌장암에서 KRAS 변이 기반 합성치사 후보를 찾고 논문으로 검증해 줘")
```

---

## 실행 예시: 왜 반증 검사가 필요한가

췌장암 KRAS 분석 실행 결과(DepMap 로컬 스냅샷, 세포주 47종 = KRAS 변이 40 / 야생형 4 / 미프로파일 3):

통계 필터를 통과한 상위 후보

| 유전자 | 변이군 | 야생형 | 차이 | q | 변이 의존 |
|---|---|---|---|---|---|
| VPS4A | −0.661 | −0.179 | −0.482 | 0.005 | 52% |
| NDE1 | −0.649 | −0.308 | −0.341 | 0.047 | 62% |
| SREBF1 | −0.547 | −0.215 | −0.332 | 0.029 | 42% |

**통계상 1순위였던 VPS4A는 `check_dependency_confounders`에서 탈락합니다.**

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
