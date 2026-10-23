# transpose.awk - swap rows and columns of a whitespace separated table
# usage: awk -f transpose.awk table.txt
{ for (i = 1; i <= NF; i++) a[i, NR] = $i; if (NF > maxnf) maxnf = NF }
END {
    for (i = 1; i <= maxnf; i++) {
        line = ""
        for (j = 1; j <= NR; j++) line = line (j > 1 ? " " : "") a[i, j]
        print line
    }
}
