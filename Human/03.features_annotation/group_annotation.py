# Step 1: 从文件3创建原始注释到组的映射
annotation_to_group = {}
with open('refData/state_annotations_processed', 'r') as f:
    # skip header
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        annotation = parts[1]
        group = parts[2]
        annotation_to_group[annotation] = group

all_groups = sorted(set(annotation_to_group.values()))

# Step 2: 读取文件1并为文件2中的每个区域进行注释，计算各组的占比
annotations = {}
with open('refData/hg38lift_genome_100_segments.bed', 'r') as f:
    for line in f:
        parts = line.strip().split()
        chr_name = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        annotation = parts[3].split('_')[-1]  # 只保留下划线之后的部分
        group = annotation_to_group.get(annotation)

        if group:  # 只有在知道组的情况下才添加注释
            if chr_name not in annotations:
                annotations[chr_name] = []
            annotations[chr_name].append((start, end, group))

group_percentages = {}
with open('query.list', 'r') as f:
    f.readline()
    for line in f:
        parts = line.strip().split(':')
        chr_name = parts[0]
        start, end = map(int, parts[1].split('-'))
        total_length = end - start
        
        group_counts = {g: 0 for g in all_groups}
        if chr_name in annotations:
            for a_start, a_end, group in annotations[chr_name]:
                if a_start < end and a_end > start:
                    overlap_start = max(start, a_start)
                    overlap_end = min(end, a_end)
                    overlap_length = overlap_end - overlap_start
                    group_counts[group] += overlap_length

        region_name = f"{chr_name}:{start}-{end}"
        percentages = {g: (group_counts[g] / total_length) * 100 for g in all_groups}
        group_percentages[region_name] = percentages

# Step 3: 将各组的占比写入新的输出文件
with open('grouped_output.txt', 'w') as out:
    # Write header
    out.write("Region\t" + "\t".join(all_groups) + "\n")
    
    for region, percentages in group_percentages.items():
        out.write(region + "\t" + "\t".join([str(percentages[g]) for g in all_groups]) + "\n")


