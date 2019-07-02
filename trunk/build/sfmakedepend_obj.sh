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
    echo "   ${script_name} src_dir obj_dir handle file"
    echo " "
    echo " "
    
    # Exit with non-zero exit status
    exit 1
}

main()
{
    # Check the number of arguments after command line option processing.
    if [ $# != "4" ]; then
	print_usage
    fi
    
    # Extract arguments and initialize.
    script_name=${0##*/}
    src_dir="$1"
    obj_dir="$2"
    handle="$3"
    file="$4"
    out_file_obj=${handle}_obj.mk
    out_file_dep=${handle}_dep.mk
    one_tab="\t"    

    # Check existence of the directory provided through the command-line.
    if [ ! -d "${src_dir}" ]; then
	echo "${script_name}: Source directory does not exist (${src_dir})."
	exit 1
    fi

    # Check existence of the file provided through the command-line
    if [ ! -f "${file}" ]; then
	echo "${script_name}: File does not exist (${file})."
	exit 1
    fi
    
    src_files=""
    # Read file line by line
    while read fline; do
        if [ -d "${src_dir}/${fline}" ]; then
            
            # Get a list of the files in the source directory.
            list_mks="$(find ${src_dir}/${fline} -name "list.mk" -print)"

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
    done < ${file}

    obj_files=""
    echo "  " > ${out_file_dep}
    for src_file in ${src_files}; do
        obj_file=${src_file##*/}
        obj_files="${obj_files} ${obj_dir}/${obj_file%.*}.o"
        case ${src_file##*.} in
            "F")
                echo "${obj_dir}/${obj_file%.*}.o : ${src_file} " >> ${out_file_dep}
                echo -e "${one_tab}@echo \"Compiling ${obj_dir}/${src_file}\"" >> ${out_file_dep}
                echo -e "${one_tab}\$(FC_WORK) -o \$@ -c \$<" >> ${out_file_dep}
                ;;
            "F90")
                echo "${obj_dir}/${obj_file%.*}.o : ${src_file} " >> ${out_file_dep}
                echo -e "${one_tab}@echo \"Compiling ${obj_dir}/${src_file}\"" >> ${out_file_dep}
                echo -e "${one_tab}\$(FC_WORK) -o \$@ -c \$<" >> ${out_file_dep}
                ;;
            "c")
                echo "${obj_dir}/${obj_file%.*}.o : ${src_file} " >> ${out_file_dep}
                echo -e "${one_tab}@echo \"Compiling ${obj_dir}/${src_file}\"" >> ${out_file_dep}
                echo -e "${one_tab}\$(CC_WORK) -o \$@ -c \$<" >> ${out_file_dep}
        esac
    done
    echo -e "${handle}_OBJ += ${obj_files}" > ${out_file_obj}
    
    return 0
}

main "$@"
