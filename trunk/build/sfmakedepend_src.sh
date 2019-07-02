#!/bin/bash
#

print_usage()
{
    local script_name
    
    # Get the script name
    script_name=${0##*/}
    
    # Echo usage info
    echo " "
    echo " "${script_name}
    echo " "
    echo " Usage:"
    echo "   ${script_name} src_dir  out_file"
    echo " "
    echo " "
    
    # Exit with non-zero exit status
    exit 1
}

main()
{
    # Check the number of arguments after command line option processing.
    if [ $# != "2" ]; then
	print_usage
    fi
    
    # Extract arguments and initialize.
    script_name=${0##*/}
    src_dir="$1"
    out_file="$2"
    sub_dir="common.dir gmp.dir hp3d.dir"
    one_tab="\t"    

    # Check existence of the directory provided through the command-line.
    if [ ! -d "${src_dir}" ]; then
	echo "${script_name}: Source directory does not exist (${src_dir})."
	exit 1
    fi

    # Check that the file provided through the command-line exists.                             
    rm -f ${out_file}
    echo "${script_name}: searching ${sub_dir}"
    for sub in ${sub_dir}; do
        dir_file=${src_dir}/${sub}
        if [ -f ${dir_file} ]; then
            echo "${script_name}: ${dir_file} is detected."
        else
            echo "${script_name}: ${dir_file} is NOT detected."
            exit 1
        fi
        src_files=""
        while read fline; do

            # Check that the subdirectory exist.
            if [ -d "${src_dir}/${fline}" ]; then
                list_mks="$(find ${src_dir}/${fline} -name "list.mk" -print)"

                # Gather all list.mk
                for list_mk in ${list_mks}; do
                    while read line; do
                        src_file=${list_mk%/*}/${line}
                        if [ -f ${src_file} ]; then
                            src_files="${src_files} ${src_file}"
                        elif [ ${line} != "0" ]; then
                            echo "${script_name}: Source file does not exist (${src_file})."
                            exit 1
                        fi
                    done < ${list_mk}
                done
            elif [ ${fline} != "0" ]; then
                echo "${script_name}: Source directory does not exist (${src_dir}/${fline})."
                exit 1
            fi
        done < ${dir_file}

        dir_name=$(echo ${sub} | sed 's/\.[^.]*$//')
        echo -e "MK_${dir_name^^}_SRC += ${src_files}" >> ${out_file}
    done    

    return 0
}

main "$@"
