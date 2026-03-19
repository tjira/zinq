#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "USAGE: $0 <INPUT> <QUEUE>" && exit 1
fi

echo -e "#!/bin/bash\n\nzinq $1" > ${1%.*}_run.sh && chmod +x ${1%.*}_run.sh && qsub -cwd -q "$2" -V ${1%.*}_run.sh
