grep -v '>' $1 | awk '!seen[$0]++'
