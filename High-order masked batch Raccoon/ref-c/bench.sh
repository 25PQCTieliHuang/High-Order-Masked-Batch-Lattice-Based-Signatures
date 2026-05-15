#!/bin/bash

MODE="batch"

while [ $# -gt 0 ]; do
    case "$1" in
        -b|--batch)
            MODE="batch"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [-b|--batch]"
            echo "  -b, --batch    Batch mode (save logs, minimal output)"
            echo "  -h, --help     Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

for dut in \
    RACCOON_128_1   RACCOON_128_2   RACCOON_128_4   \
    RACCOON_128_8   RACCOON_128_16  RACCOON_128_32  \
    RACCOON_192_1   RACCOON_192_2   RACCOON_192_4   \
    RACCOON_192_8   RACCOON_192_16  RACCOON_192_32  \
    RACCOON_256_1   RACCOON_256_2   RACCOON_256_4   \
    RACCOON_256_8   RACCOON_256_16  RACCOON_256_32
do
    logf="bench_${dut}.txt"
    echo "=== $dut [MODE: $MODE] ==="

    # 清理并编译
    make obj-clean
    if [ "$MODE" = "batch" ]; then
        make MODE="$MODE" RACCF="-D${dut} -DBENCH_TIMEOUT=2.0"
    else
        make RACCF="-D"$dut" -DBENCH_TIMEOUT=2.0"
    fi

    if [ ! -x ./xtest ]; then
        echo "ERROR: xtest not executable" | tee -a "$logf"
        continue
    fi

    if [ "$MODE" = "batch" ]; then
        # batch 模式：输出到文件，终端只显示简要信息
        {
            echo "=== $logf ==="
            
            echo "MODE: $MODE"
            echo "---"
        } > "$logf"                # 先写入文件头

        # 将 xtest 的输出同时显示在终端？batch 模式可以隐藏，只写入文件
        ./xtest 2>&1 | tee -a "$logf" > /dev/null 2>&1
        # 或者如果你完全不想在终端看到输出，可以用：
        # ./xtest >> "$logf" 2>&1

        {
            echo "---"
            echo "STACK_USAGE:"
            if [ -f racc_core.su ]; then
                grep -E '_core_keygen|_core_sign|_core_verify' racc_core.su
            else
                echo "racc_core.su not found"
            fi
        } >> "$logf" 2>&1

        echo "[Saved to $logf]"
    else
        
        {
            echo "=== $logf ==="
            echo "MODE: $MODE"
            echo "---"
            ./xtest 2>&1
            echo "---"
            echo "STACK_USAGE:"
            if [ -f racc_core.su ]; then
                grep -E '_core_keygen|_core_sign|_core_verify' racc_core.su
            else
                echo "racc_core.su not found"
            fi
        } | tee "$logf"           
    fi
done