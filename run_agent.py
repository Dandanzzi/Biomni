import pandas as pd
from scipy.stats import ttest_ind

# ==========================================
# 1. 합성치사 분석 툴(도구) 정의
# ==========================================
def find_synthetic_lethality_candidates(cancer_type: str, target_mutation: str) -> str:
    """
    특정 암종(예: Pancreatic Cancer)에서 특정 유전자 돌연변이(예: KRAS) 유무에 따라,
    생존에 치명적인 영향을 받는 합성치사(Synthetic Lethality) 후보 유전자를 T-test로 찾습니다.
    """
    print(f"🔬 [에이전트 내부 실행 중] {cancer_type} 세포주에서 {target_mutation} 변이 기반 합성치사 분석 시작...\n")
    
    try:
        # [데이터 준비 단계: 가상 샘플 데이터]
        cell_lines_data = {
            'Cell_Line': ['MIA-PaCa-2', 'PANC-1', 'BxPC-3', 'AsPC-1', 'Capan-2', 'HPAC'],
            'Cancer_Type': ['Pancreatic Cancer'] * 6,
            f'{target_mutation}_Mutated': [True, True, False, True, False, False]
        }
        df_cell = pd.DataFrame(cell_lines_data)
        
        dependency_data = {
            'Cell_Line': ['MIA-PaCa-2', 'PANC-1', 'BxPC-3', 'AsPC-1', 'Capan-2', 'HPAC'],
            'Gene_A_Score': [-0.1, -0.2, -0.1, -0.3, -0.2, -0.1], 
            'Gene_B_Score': [-1.5, -1.8, -0.1, -1.6, -0.2, -0.1], # KRAS 돌연변이 시 치명적
            'Gene_C_Score': [-0.9, -0.8, -1.0, -0.7, -0.9, -0.8]  
        }
        df_dep = pd.DataFrame(dependency_data)
        
        # 데이터 병합
        df_merged = pd.merge(df_cell, df_dep, on='Cell_Line')
        
        # [분석 단계: T-test 통계 검증]
        mutated_group = df_merged[df_merged[f'{target_mutation}_Mutated'] == True]
        wt_group = df_merged[df_merged[f'{target_mutation}_Mutated'] == False]
        
        results = []
        genes_to_test = ['Gene_A_Score', 'Gene_B_Score', 'Gene_C_Score']
        
        for gene in genes_to_test:
            stat, p_value = ttest_ind(mutated_group[gene], wt_group[gene])
            
            mut_mean = mutated_group[gene].mean()
            wt_mean = wt_group[gene].mean()
            
            if mut_mean < wt_mean and p_value < 0.05:
                results.append(
                    f"- 타겟 유전자: {gene.replace('_Score', '')}\n"
                    f"  * {target_mutation} 변이 세포주 평균 점수: {mut_mean:.2f} (치명적)\n"
                    f"  * {target_mutation} 정상 세포주 평균 점수: {wt_mean:.2f} (안전함)\n"
                    f"  * p-value: {p_value:.4f} (통계적 유의성 확보)"
                )
        
        # [결과 반환]
        if results:
            return f"✅ [{cancer_type}] {target_mutation} 변이 기반 합성치사 분석 결과:\n" + "\n".join(results)
        else:
            return f"[{cancer_type}]에서 {target_mutation} 변이와 유의미한 합성치사 관계를 가진 유전자를 찾지 못했습니다."
            
    except Exception as e:
        return f"데이터 분석 중 오류가 발생했습니다: {str(e)}"

# ==========================================
# 2. 에이전트 셋팅 (임시)
# (실제 Biomni 라이브러리를 사용 중이라면 이 부분을 원래 코드로 교체하세요)
# ==========================================
class MockBiomniAgent:
    def __init__(self):
        self.tools = []
        
    def add_tool(self, tool_function):
        self.tools.append(tool_function)
        print(f"⚙️  툴 등록 완료: {tool_function.__name__}")
        
    def run(self, user_query: str):
        # 챗봇이 질문을 이해하고 방금 만든 툴을 꺼내 쓰는 과정을 시뮬레이션
        if "췌장암" in user_query and "KRAS" in user_query:
            # 에이전트가 툴을 실행하고 그 결과를 반환함
            return find_synthetic_lethality_candidates(cancer_type="Pancreatic Cancer", target_mutation="KRAS")
        else:
            return "질문에 맞는 분석 툴을 찾지 못했습니다."

# 에이전트 생성 및 툴 장착
agent = MockBiomniAgent()
agent.add_tool(find_synthetic_lethality_candidates)


# ==========================================
# 3. 실행 스위치 (이 부분이 있어야 터미널에서 작동합니다!)
# ==========================================
if __name__ == "__main__":
    print("\n" + "="*50)
    
    question = "췌장암에서 KRAS 돌연변이와 합성치사 관계인 후보 유전자를 찾아줘"
    print(f"🧑‍🔬 사용자 질문: {question}\n")
    print("🤖 에이전트가 데이터 분석을 시작합니다... (잠시만 기다려주세요)\n")
    
    # 에이전트에게 질문을 던지고 답변을 받아옴
    answer = agent.run(question)
    
    print("================== [에이전트 분석 결과] ==================")
    print(answer)
    print("=========================================================\n")