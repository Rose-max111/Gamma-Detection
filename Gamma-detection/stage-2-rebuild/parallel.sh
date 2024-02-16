#!/bin/bash
mkdir -p output

# 运行任务的函数
run_task() {
    local script="$1"
    ./"$script" &
}

# 脚本数量
script_count="$1"

# 并行运行脚本
for ((i=1; i<=script_count; i++)); do
    script="script/run${i}.sh"
    run_task "$script"
done

# 等待所有后台任务完成
wait
