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
    echo "   ${script_name} src_dir obj_dir out_file_obj out_file_src"
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
    out_file_obj="$3"
    out_file_src="$4"
    one_tab="\t"    

    # Check existence of the directory provided through the command-line.
    if [ ! -d "${src_dir}" ]; then
	echo "${script_name}: Source directory does not exist (${src_dir})."
	exit 1
    fi
    
    # Get a list of the files in the source directory.
    list_mks="$(find ${src_dir} -name "list.mk" -print)"


    src_files=""
    for list_mk in ${list_mks}; do
        echo "Collecting object files from ${list_mk}"
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

    obj_files=""
    for src_file in ${src_files}; do
        case ${src_file##*.} in
            "F")
                obj_files="${obj_files} ${obj_dir}/${src_file%.F}.o" ;;
            "F90")
                obj_files="${obj_files} ${obj_dir}/${src_file%.F90}.o" ;;
            "c")
                obj_files="${obj_files} ${obj_dir}/${src_file%.c}.o" ;;
        esac 
    done
    echo -e "OBJ= ${obj_files}" > ${out_file_obj}
    echo -e "SRC_FILE = ${src_files}" > ${out_file_src}
    return 0
}

main "$@"
