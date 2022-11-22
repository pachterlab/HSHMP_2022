grep -v '>' $1 | awk '!seen[$0]++' | wc -l
