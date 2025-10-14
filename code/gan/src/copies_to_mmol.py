import pandas as pd

# the proteomic data file
file_path = "data/E.coli/proteomics.csv"
df = pd.read_csv(file_path)

# the column of "Molecular weight (Da)"
mw_col = "Molecular weight (Da)"
abundance_all = df.columns[4:]  # protein abundance( fg/cell)

# dry cell weight
dcw_per_cell = 1e-13  #  g/cell

#  Protein Mass (fg/cell) to mmol/g DCW
for col in abundance_all:
    df[col] = (df[col] * 1e-15 * 1000) / (df[mw_col] * dcw_per_cell)

# 保存为新CSV文件
output_file = "data/E.coli/proteomic_all.csv"
df.to_csv(output_file, index=False)

print(f"Finished writing {output_file}")
