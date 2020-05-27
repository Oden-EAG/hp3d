#/bin/bash
if [ $# -eq 0 ];
then
 echo "Syntax: $(basename $0) ext_from ext_to "
 exit 1
fi

find . -type f -name '*.'$1 -print | while read file; do
  echo "renaming $file to $(basename $file .$1).$2";
  mv "$file" "$(dirname $file)/$(basename $file .$1).$2";
done


