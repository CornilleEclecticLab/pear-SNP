while IFS= read -r line; do chr=$(echo "$line" | awk '{print $1}') ; echo "$line" > "${chr}.fasta.fai" ;done < "fasta.fai"
