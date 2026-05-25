import os
import shutil

samples_file = r"C:\Users\Asus\Desktop\EHI_laterality\Part4_expression\RNA_samples.txt"
source_dir = r"D:\imputed_expression_250526\Brain_Substantia_nigra\imputed_counts"
target_dir = r"C:\Users\Asus\Desktop\EHI_laterality\Part4_expression\Imputed_brain\Brain_Substantia_nigra"

with open(samples_file, "r") as f:
    valid_samples = set(
        line.strip().replace(".tsv", "")
        for line in f
        if line.strip()
    )

os.makedirs(target_dir, exist_ok=True)

copied = 0
skipped = 0

for filename in os.listdir(source_dir):

    if filename.endswith(".tsv"):
        parts = filename.split("_")
        sample_name = f"{parts[0]}_{parts[1]}"

        if sample_name in valid_samples:
            shutil.copy2(
                os.path.join(source_dir, filename),
                os.path.join(target_dir, filename)
            )
            copied += 1
            print(f"[COPIADO] {filename}")
        else:
            skipped += 1

print("\nProceso terminado")
print(f"Archivos copiados: {copied}")
print(f"Archivos ignorados: {skipped}")