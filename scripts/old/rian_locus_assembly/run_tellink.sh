#!/bin/bash

if [ "$#" == "0" ]; then
echo "Error: No arguments" >&2; exit 1;
fi

if [ "$#" == "1" ]; then
        if [ "$1" == "-v" ] || [ "$1" == "--version" ]; then
                docker run docker-tellink cat /tellysis/tellink/version.txt; exit 1;
        fi
fi

while (( "$#" )); do
  case "$1" in
    -r1|--read1)
      r1=$2
      shift 2
      ;;
    -r2|--read2)
      r2=$2
      shift 2
      ;;
    -i1|--i1)
      i1=$2
      shift 2
      ;;
    -o|--output)
      output_path=$2
      shift 2
      ;;
    -r|--reference)
      r=$2
      shift 2
      ;;
    -d|--metagenomics)
      d=$2
      shift 2
      ;;
    -p|--prefix)
      pref=$2
      shift 2
      ;;
    -k|--kmer)
      kmer=$2
      shift 2
      ;;
    -lc|--local-kmer)
      lc=$2
      shift 2
      ;;
    -j|--cores)
      cores=$2
      shift 2
      ;;
    -*|--*=|*) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

mkdir -p $output_path
cd $output_path

#get host system's directories for volume mapping
i1_path=""
r1_path=""
r2_vol=""
r2_opt=""
r_vol=""
r_opt=""
d_opt=""

i1_path=`readlink -f $(dirname $i1)`
r1_path=`readlink -f $(dirname $r1)`
if [ "$r2" != '' ]; then
	r2_path=`readlink -f $(dirname $r2)`
	r2_vol="-v $r2_path:$r2_path"
	r2_opt="-r2 $r2"
fi

if [ "$r" != '' ]; then
	r_path=`readlink -f $(dirname $r)`
	r_vol="-v $r_path:$r_path"
	r_opt="-r $r"
fi

if [ "$d" != '' ]; then
        d_opt="-d $d"
fi


#get host user for mapping
uid=$(id -u ${USER})
gid=$(id -g ${USER})

#get output directory for mapping
outp=`readlink -f $output_path`
mkdir -p $outp

echo "$USER:x:$uid:$gid::/home/ubuntu:/bin/bash" > $outp/temp_1
echo "$USER:x:$gid:" > $outp/temp_2

#run tellink docker image
docker run \
	-v $i1_path:$i1_path \
	-v $r1_path:$r1_path $r2_vol \
	$r_vol \
	-v $outp:$output_path \
	--mount type=bind,source=$output_path/temp_1,target=/etc/passwd \
	--mount type=bind,source=$output_path/temp_2,target=/etc/group \
	-u $uid:$gid  \
	--rm \
        018daaacd197 \
	bash -c "/tellysis/tellink/StartTellinkPipeline -r1 $r1 $r2_opt -i1 $i1 $r_opt $d_opt \
            -o $output_path -p $pref -k $kmer -lc $lc -l ust -j $cores \
            2>&1 | tee $output_path/$(basename $output_path).log"

# a737b0c2fbc4 (Docker image with modifications)
# 0fe94f4429fa (Docker image with the use of an extension of QUAST for large-scale genome assemblies)
# 85b0e8a6477c (Original Docker image)