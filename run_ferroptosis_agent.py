"""바이옴니 에이전트에게 'SREBF1 결손이 췌장암에서 페롭토시스를 일으키는가'를 직접 조사시킨다.

워크플로를 지정하지 않습니다. 질문만 던지면 Biomni의 tool retriever가 필요한 도구를 고르고,
에이전트가 <execute> 파이썬 블록 안에서 DepMap 분석과 PubMed 조회를 스스로 수행합니다.

    python run_ferroptosis_agent.py                       # 기본 질문
    python run_ferroptosis_agent.py --query "다른 질문"    # 임의 질문
    python run_ferroptosis_agent.py --save answer.md      # 최종 답변 저장

LLM API 키가 필요합니다 (ANTHROPIC_API_KEY 등). 실행마다 과금됩니다.
"""

import argparse
import os

from dotenv import load_dotenv

DEFAULT_QUESTION = """췌장암(PDAC)에서 SREBF1 결손이 세포를 죽이는 기전이 페롭토시스(ferroptosis)인지 조사해 줘.

확인하고 싶은 것:
1. DepMap CRISPR 데이터에서 SREBF1 의존성이 페롭토시스 관련 유전자(SCD, GPX4, ACSL4, SLC7A11,
   AIFM2, LPCAT3)의 의존성과 어떻게 연관되는가. SREBP 경로 자체(SCAP, MBTPS1, MBTPS2)와의
   상관은 분석이 제대로 작동하는지 확인하는 양성 대조로 함께 봐 줘.
2. KRAS 변이 췌장암 세포주에서 이 유전자들의 실제 의존도는 얼마인가.
3. PubMed에 SREBF1-페롭토시스-췌장암을 연결하는 근거가 있는가. 반대로 이를 부정하는 보고는 없는가.
4. 정상 조직에서의 발현을 근거로 치료 window가 있다고 볼 수 있는가.
5. 위 근거를 종합해, 페롭토시스 억제제(ferrostatin-1)로 구제되는지를 검증할 실험을 설계해 줘.

중요: 근거가 약한 부분은 약하다고 명시하고, 지지 근거뿐 아니라 반대 근거와 교란 가능성도
함께 보고해 줘. 확정되지 않은 것을 확정된 것처럼 서술하지 마."""


def main() -> None:
    parser = argparse.ArgumentParser(description="Biomni 에이전트에게 SREBF1-페롭토시스 가설을 조사시킨다")
    parser.add_argument("--query", default=DEFAULT_QUESTION, help="에이전트에게 던질 질문")
    parser.add_argument("--llm", default="claude-sonnet-4-5-20250929")
    parser.add_argument("--data-path", default="./data")
    parser.add_argument("--save", default=None, help="최종 답변을 저장할 파일 경로")
    args = parser.parse_args()

    load_dotenv(".env")
    if not any(os.getenv(k) for k in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "AZURE_OPENAI_API_KEY")):
        raise SystemExit("LLM API 키가 없습니다. ANTHROPIC_API_KEY 등을 먼저 설정하세요.")

    from biomni.agent import A1

    agent = A1(path=args.data_path, llm=args.llm)
    registered = [t["name"] for t in agent.module2api.get("biomni.tool.synthetic_lethality", [])]
    registered += [t["name"] for t in agent.module2api.get("biomni.tool.organoid_sl", [])]
    print(f"등록된 합성치사/오가노이드 도구: {registered}")
    print("도구 선택과 실행 순서는 에이전트가 스스로 계획합니다.\n")
    print("=" * 78)
    print("질문")
    print("=" * 78)
    print(args.query)
    print("")

    _, answer = agent.go(args.query)

    print("\n" + "=" * 78)
    print("최종 답변")
    print("=" * 78)
    print(answer)

    if args.save:
        with open(args.save, "w") as handle:
            handle.write(answer)
        print(f"\n최종 답변을 {args.save}에 저장했습니다.")


if __name__ == "__main__":
    main()
