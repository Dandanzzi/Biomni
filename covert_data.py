import pandas as pd

# TSV 파일 읽기
file_path = "/data2/home/delee/biomni/data/biomni_data/data_lake/gene_sl_gene.tsv"
raw = pd.read_csv(file_path, sep='\t')

# 확인된 컬럼명(x_name, y_name)을 사용하여 정확하게 매핑
table = pd.DataFrame({
    "gene_a": raw["x_name"].astype(str).str.upper().str.strip(),
    "gene_b": raw["y_name"].astype(str).str.upper().str.strip(),
    "score": pd.to_numeric(raw["r.statistic_score"], errors="coerce").fillna(1.0) if "r.statistic_score" in raw.columns else 1.0,
    "source": raw["rel_source"].astype(str) if "rel_source" in raw.columns else "",
})

# 동일 유전자 쌍 제거 및 결측치(문자열 "NAN" 포함) 제거
table = table[(table.gene_a != table.gene_b) & (table.gene_a != "NAN") & (table.gene_b != "NAN")]

# Parquet 파일로 저장
output_path = "/data2/home/delee/biomni/data/biomni_data/data_lake/synlethdb_human_sl.parquet"
table.to_parquet(output_path, index=False)
print(f"변환 완료! {len(table)}개의 유전자 쌍이 성공적으로 저장되었습니다.")

print(table.head())