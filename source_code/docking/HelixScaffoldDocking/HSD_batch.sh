#!/bin/bash

# 解析命令行参数
while getopts "t:b:a:" opt; do
  case $opt in
    t) target="$OPTARG"
    ;;
    b) batch_info="$OPTARG"
    ;;
    a) atomid="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# 检查是否所有参数都已提供
if [ -z "$target" ] || [ -z "$batch_info" ] || [ -z "$atomid" ]; then
  echo "Usage: $0 -t target -b batch_info -a atomid"
  exit 1
fi

# 执行命令
$HSD/docking/HelixScaffoldDocking/HelixScaffoldDocking_batch "$target" "$batch_info" "$atomid"