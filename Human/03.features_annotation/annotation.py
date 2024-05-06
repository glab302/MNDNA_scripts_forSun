# 读取文件1，建立注释字典
annotations = {}
all_annotation_types = set()

with open('refData/hg38lift_genome_100_segments.bed', 'r') as f:
    for line in f:
        parts = line.strip().split()
        chr_name = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        annotation = parts[3]
        
        all_annotation_types.add(annotation)
        
        if chr_name not in annotations:
            annotations[chr_name] = []
        annotations[chr_name].append((start, end, annotation))

# 根据前缀进行排序
all_annotation_types = sorted(list(all_annotation_types), key=lambda x: (int(x.split('_')[0]), x))

# 创建一个输出文件
with open('output.txt', 'w') as out:
    # 写入表头
    out.write("Region\t" + "\t".join(all_annotation_types) + "\n")
    
    # 处理文件2中的区域，查找其注释并计算占比
    with open('query.list', 'r') as f:
        # skip header
        f.readline()
        for line in f:
            parts = line.strip().split(':')
            chr_name = parts[0]
            start, end = map(int, parts[1].split('-'))
            total_length = end - start

            annotation_counts = {atype: 0 for atype in all_annotation_types}
            if chr_name in annotations:
                for a_start, a_end, annotation in annotations[chr_name]:
                    if a_start < end and a_end > start:
                        overlap_start = max(start, a_start)
                        overlap_end = min(end, a_end)
                        overlap_length = overlap_end - overlap_start
                        annotation_counts[annotation] += overlap_length

            percentages = [str((annotation_counts[atype] / total_length) * 100) for atype in all_annotation_types]
            out.write(f"{chr_name}:{start}-{end}\t" + "\t".join(percentages) + "\n")


