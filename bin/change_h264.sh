#!/bin/bash

if [[ $1 -eq 1 ]]
then
file_name=$2
else
file_name="bus_cavlc_multi_idr"
fi

echo $file_name

raw_h264_file=${file_name}"_Copy.264"
new_h264_file=${file_name}".264"

v_dir="../vediofile/decoder"

rm -f ${v_dir}"/"${new_h264_file}
cp ${v_dir}"/"${raw_h264_file} ${v_dir}"/"${new_h264_file}

echo "cp "${raw_h264_file}" to "${new_h264_file}
