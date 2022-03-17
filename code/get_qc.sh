
#!/bin/sh

## use this script to create a spreadsheet with file name and counts (either numb peaks or numb reads)
## from Davide Vespasiani

helpFunction()
{
   echo ""
   echo "Usage: $0 -g genome -i input_dir -f file_type"
   echo -e "\t-g genome assembly (e.g. hg38)"
   echo -e "\t-i Directory where the bams/narrowPeak files are stored"
   echo -e "\t-f file type that must be processed (e.g., bam)"
   exit 1 # Exit script after printing help
}

while getopts "g:i:f:" flag; do
    case "${flag}" in
        g) genome="$OPTARG";;
        i) input_dir="$OPTARG";;
        f) file_type="$OPTARG";;
        ?) helpFunction ;; # Print helpFunction 
    esac
done


basedir='/data/projects/punim0586/lecook/chipseq-cross-species/'
wd="$basedir/output/$input_dir"

echo "Working dir is $wd"

file_dir="$wd/bam_files"
qc_dir="$wd/qc"
out_file="$qc_dir/finalqc.txt"

## run command
echo "Getting filenames and count"

for f in "$file_dir/"*."${file_type}"; do
        filname=${f#"$file_dir/"}
        printf '%s\n' "$filname"  >> names.out
    if [[ $f == *.gz ]]; then
        gunzip -nc "$f" | wc -l >> counts.out
    elif [[ $f == *.bam ]]; then
        module load samtools/1.9
        samtools view "$f" | wc -l >> counts.out 
    fi
done

cd "$wd/Files"

paste -d "\t" names.out counts.out > "$qc_dir"/stat_qc.out
rm names.out 
rm counts.out 

## Call custom python script to make final spreadsheet

echo "Making spreadsheet"
python "$basedir"/code/make_spreadsheet.py "$qc_dir"/stat_qc.out $out_file

rm "$qc_dir"/stat_qc.out