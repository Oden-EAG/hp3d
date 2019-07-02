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
    echo "   ${script_name} src_dir obj_dir mod_dir out_file"
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
    mod_dir="$3"
    out_file="$4"
    one_tab="\t"    

    # If out_file exist, remove
    if [ -f "${out_file}" ]; then
        rm -f ${out_file}
    fi

    # Check existence of the directory provided through the command-line.
    if [ ! -d "${src_dir}" ]; then
	echo "${script_name}: Source directory does not exist (${src_dir})."
	exit 1
    fi
    
    # Get a list of the files in the source directory.
    list_mks="$(find ${src_dir} -name "list.mk" -print)"

    src_files=""
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

    # Extract modules
    mod_files=""
    for src_file in ${src_files}; do
        mod=${src_file##*/}
        mod=${mod%.*}
        mod_check=""
        mod_check=$(grep -w -h module[[:space:]]${mod} ${src_file})
        if [ "${mod_check}" != "" ]; then
            mod_files="${mod_files} ${src_file}"
        fi
    done

    # Check whether src_files use module
    for src_file in ${src_files}; do
        modules=$(grep -w -h use ${src_file})
        prev=""
        routines=""
        for module in ${modules}; do
            if [ "${prev:0:3}" == "use" ]; then
                count=1
                for routine in ${routines}; do
                    if [ "${module}" == "${routine}" ]; then
                        count=0
                        break
                    fi
                done
                if [ ${count} == 1 ]; then
                    routines="${routines} ${module}"
                fi
            fi
            prev=${module}
        done
        
        dep_files=""
        for routine in ${routines}; do 
            count=0
            for dep_file in ${mod_files}; do
                dep=${dep_file##*/}
                dep=${dep%.*}
                if [ "${routine}" == "${dep}" ]; then
                    count=1
                    break
                fi
            done
            if [ ${count} == 1 ]; then
                dep_files="${dep_files} ${mod_dir}/${dep}.mod"
            fi
        done

        echo "Creating dependency for ${src_file}"

        case ${src_file##*.} in
            "F")
                echo "${obj_dir}/${src_file%.F}.o : ${src_file} ${dep_files}" >> ${out_file} 
                echo -e "${one_tab}@mkdir -p \$(dir \$@)" >> ${out_file}
                echo -e "${one_tab}@echo \"Compiling ${obj_dir}/${src_file}\"" >> ${out_file}
                echo -e "${one_tab}\$(FC_WORK) -o \$@ -c \$<" >> ${out_file}
                ;;
            "F90")
                echo "${obj_dir}/${src_file%.F90}.o : ${src_file} ${dep_files}" >> ${out_file} 
                echo -e "${one_tab}@mkdir -p \$(dir \$@)" >> ${out_file}
                echo -e "${one_tab}@echo \"Compiling ${obj_dir}/${src_file}\"" >> ${out_file}
                echo -e "${one_tab}\$(FC_WORK) -o \$@ -c \$<" >> ${out_file}
                ;;
            "c")
                echo "${obj_dir}/${src_file%.c}.o : ${src_file} ${dep_files}" >> ${out_file} 
                echo -e "${one_tab}@mkdir -p \$(dir \$@)" >> ${out_file}
                echo -e "${one_tab}@echo \"Compiling ${obj_dir}/${src_file}\"" >> ${out_file}
                echo -e "${one_tab}\$(CC_WORK) -o \$@ -c \$<" >> ${out_file}
        esac 
        echo -e "${one_tab}" >> ${out_file}
    done
    return 0
}

main "$@"
