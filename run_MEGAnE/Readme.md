
Usage of s02.*.sh

```bash
for path in ../output/s01.*/*; do
    popbase=$(basename "$path")
    realpath "$path"/* > "../input/s02.pop_path.$popbase.txt"
    sbatch s02.merge_no-reference_ME_insertions_and_reference_ME_polymorphisms.sh "../input/s02.pop_path.$popbase.txt" "$popbase" 
done
