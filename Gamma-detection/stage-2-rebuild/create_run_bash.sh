#!/bin/bash

# 设置总数据范围和每个脚本测试的数据数量
total_data=2000
data_per_script=$1
script_count=$((total_data / $1))

# 创建脚本目录
mkdir -p script
for ((i = 1; i <= script_count; i++)); do
    start_num=$(( (i - 1) * data_per_script ))
    end_num=$(( start_num + data_per_script - 1 ))

    script_file="script/run${i}.sh"

    cat > "$script_file" <<EOF
#!/bin/bash

start_num=$start_num
end_num=$end_num


# 循环枚举数字
for ((x = start_num; x <= end_num; x++)); do
    echo "Processing file: \${x}.txt"

    if [[ "\$x" == 671 || "\$x" == 741 || "\$x" == 1472 || "\$x" == 1647 || "\$x" == 1984 ]]; then
        echo "Skipping file: \${x}.txt"
        touch "output/\${x}_ans.txt"
    elif [ ! -f "output/\${x}_ans.txt" ]||[ ! -s "output/\${x}_ans.txt" ]; then
	# 运行 rank，并将输入输出重定向到对应的文件
        if [[ "\$x" == 142 || "\$x" == 184 || "\$x" == 228 || "\$x" == 241 || "\$x" == 277 || "\$x" == 804 || "\$x" == 902 || "\$x" == 962 || "\$x" == 1268 || "\$x" == 1567 || "\$x" == 1615 || "\$x" == 1738 ]]; then
            echo "Processing file with rank_90: \${x}.txt"
            ./rank_90 < "data/\${x}.txt" > "output/\${x}_ans.txt"
        elif [[ "\$x" == 1433 ]]; then
            echo "Processing file with rank_140: \${x}.txt"
            ./rank_140 < "data/\${x}.txt" > "output/\${x}_ans.txt"
        else
            echo "Processing file with rank_150: \${x}.txt"
            ./rank_150 < "data/\${x}.txt" > "output/\${x}_ans.txt"
        fi
    fi
done
EOF

    chmod +x "$script_file"
done
