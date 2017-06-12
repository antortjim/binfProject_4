cp $2.tsv $2_flipped.tsv
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "Text read from file: $line"
    # extract the record as a string
    RECORD=$(sed $line'q;d' $2_flipped.tsv)

    # flip the string
    FLIPPED=$(echo $RECORD | tr A "|" | tr C "@" | tr G "#" | tr T "~" | \
                   tr "|" T | tr "@" G  | tr "#" C | tr "~" A)

    # replace the old string with the new one on the spot
    sed -i $line"s/.*/$FLIPPED/" $2_flipped.tsv

    cat $2_flipped.tsv | tr " " "\t" > temp2
    mv temp2 $2_flipped.tsv

done < "$1"


