
#bash script to gleaning the data

OUTPUT_PATH=$1
PREFIX=$2
DATA=$3

zcat $DATA | awk '(NR>=7)' |cut -f 1,2 > $OUTPUT_PATH/snps.txt

(echo "LG";seq 1 30) > $OUTPUT_PATH/${PREFIX}_sepchr_lod.txt

# make a list of the numberof maps and loop over it 

ls ${OUTPUT_PATH}/${PREFIX}.*.map | while read file; do
  lod=$(basename "$file" | cut -d. -f2)
  echo "Processing LOD ${lod} from file ${file}"
  if [ ! -f "${file}" ]; then
    echo "  Warning: File ${file} not found, skipping."
    continue
  fi
  # if file is only containing 0 then skip it
  if cat "${file}" | awk 'NR>1 && $1>0{exit 1}' > /dev/null; then
    echo "  Warning: File ${file} contains only 0, skipping."
    continue
  fi
  paste "${file}" "${OUTPUT_PATH}/snps.txt" | awk '($1>0){print $1, $2}' | sort -n | uniq -c | sort -rn | head -n30 | awk -v lod="${lod}" 'BEGIN{print "lodLimit_"lod};{print $0}' | paste "${OUTPUT_PATH}/${PREFIX}_sepchr_lod.txt" - > "${OUTPUT_PATH}/${PREFIX}_sepchr_lod_prel.txt" && mv "${OUTPUT_PATH}/${PREFIX}_sepchr_lod_prel.txt" "${OUTPUT_PATH}/${PREFIX}_sepchr_lod.txt"
done

echo "Finished"
date

#check resultfile in r

