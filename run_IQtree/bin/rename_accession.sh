while read key value; do
    sed -i "s/\b$key\b/$value/g" A
done < B
