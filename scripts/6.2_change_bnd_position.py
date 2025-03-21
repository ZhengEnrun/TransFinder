import argparse

# 定义命令行参数解析
def parse_args():
    parser = argparse.ArgumentParser(description="Adjust and combine translocation data.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file.")
    return parser.parse_args()

# 主函数
def process_file(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue  # 跳过空行

            # 分割行内容
            chr1, chr2, orientation, pos1, pos2, event = line.split("\t")

            # 调整格式
            if orientation == "+-":
                new_orientation = "-+"
            elif orientation == "-+":
                new_orientation = "+-"
            else:
                new_orientation = orientation  # 保持 "++" 或 "--" 不变

            # 交换 chr 和位置
            new_line = f"{chr2}\t{chr1}\t{new_orientation}\t{pos2}\t{pos1}\t{event}\n"

            # 写入原始行和调整后的行
            outfile.write(line + "\n")  # 原始行
            outfile.write(new_line)     # 调整后的行

    print(f"文件处理完成，结果已保存到 {output_file}")

# 执行脚本
if __name__ == "__main__":
    args = parse_args()
    process_file(args.input, args.output)
