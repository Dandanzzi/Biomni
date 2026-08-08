from dotenv import load_dotenv

load_dotenv(".env")
from biomni.agent import A1  # noqa: E402

agent = A1(path="./data", llm="claude-sonnet-4-5-20250929", expected_data_lake_files=[])
log, answer = agent.go("췌장암과 합성치사 관계에 있는 것을 알려줘.")
print("\n\n########## FINAL ANSWER ##########\n")
print(answer)
