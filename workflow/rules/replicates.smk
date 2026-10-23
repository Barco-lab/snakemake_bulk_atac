# 重复一致性：组内 peak 重现性 -> consensus peak 计数 -> 样本相关性、PCA；性别检查（chrX / chrY reads 占比）
# 分组见 samples.tsv 的 group 列（为空时按样本名去掉 _repN 推断）


def replicates_on():
    return REPLICATES["enabled"] and len(SAMPLES) >= 2


def replicate_figures():
    """报告中的图（不含扩展名）；条件不满足的图不画。"""
    f = ["correlation"]
    if len(SAMPLES) >= 3:
        f.append("pca")
    if REP_GROUPS:
        f.append("peak_overlap")
    if venn_groups():
        f.append("peak_venn")
    if REPLICATES["sex_check"]:
        f.append("sex")
    return f


def venn_groups():
    """画 Venn 图的组：2–4 个重复（更多时 Venn 难以辨认，只看 Jaccard）。"""
    return [g for g in REP_GROUPS if len(GROUPS[g]) <= 4]


def replicate_figure_files():
    return expand(
        "results/qc/replicates/figures/{fig}.{ext}", fig=replicate_figures(), ext=["png", "pdf", "svg"]
    )


def replicate_report_inputs():
    """QC 报告「重复一致性」一节的输入；未启用时为空。"""
    if not replicates_on():
        return {}
    d = "results/qc/replicates"
    return dict(
        rep_info=f"{d}/info.json",
        rep_samples=f"{d}/samples.tsv",
        rep_groups=f"{d}/groups.tsv",
        rep_svg=expand(f"{d}/figures/{{fig}}.svg", fig=replicate_figures()),
    )


def region_len(region):
    """chr:start-end（1-based，含两端）的长度。"""
    start, end = region.split(":")[1].replace(",", "").split("-")
    return int(end) - int(start) + 1


rule sex_counts:
    """全部比对 reads、chrX 与 chrY 计数区（references 中的 y_region）上的 reads 数（过滤去重后的 BAM）。"""
    input:
        bam="results/02_align/{sample}.dedup.bam",
        bai="results/02_align/{sample}.dedup.bam.bai",
    output:
        "results/qc/replicates/sex/{sample}.tsv",
    params:
        x=REF["chrX"],
        y=REF["y_region"],
        y_len=region_len(REF["y_region"]),
    conda:
        ENV_MACS
    shell:
        """
        y=$(samtools view -c {input.bam} {params.y})
        samtools idxstats {input.bam} | awk -v OFS='\\t' -v s={wildcards.sample} -v x={params.x} -v y=$y -v yl={params.y_len} '
            {{ t += $3 }}
            $1 == x {{ xr = $3 }}
            END {{
                print "sample", "total_reads", "x_reads", "y_reads", "y_len"
                print s, t, xr + 0, y, yl
            }}' > {output}
        """


rule reproducible_peaks:
    """组内各重复的 peak 合并（bedtools merge），保留至少在 min_reps 个重复中出现的区域；第 4 列为重复数。"""
    input:
        lambda wc: [peak_file(s) for s in GROUPS[wc.group]],
    output:
        "results/03_peaks/reproducible/{group}.bed",
    params:
        names=lambda wc: " ".join(GROUPS[wc.group]),
        n=lambda wc: min(REPLICATES["min_reps"], len(GROUPS[wc.group])),
    conda:
        ENV_MACS
    shell:
        """
        names=({params.names}); files=({input})
        for i in "${{!files[@]}}"; do
            awk -v OFS='\\t' -v n="${{names[$i]}}" '{{ print $1, $2, $3, n }}' "${{files[$i]}}"
        done | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o count_distinct \
            | awk -v n={params.n} '$4 >= n' > {output}
        """


rule consensus_peaks:
    """各组 reproducible peak 的并集；只有一个样本的组用该样本的 peak。"""
    input:
        rep=expand("results/03_peaks/reproducible/{group}.bed", group=REP_GROUPS),
        single=[peak_file(GROUPS[g][0]) for g in GROUPS if g not in REP_GROUPS],
    output:
        "results/03_peaks/consensus.bed",
    conda:
        ENV_MACS
    shell:
        """
        cat {input.rep} {input.single} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > {output}
        """


rule peak_counts:
    """每个样本在 consensus peak 上的计数（featureCounts；PE 按片段计）。"""
    input:
        bam="results/02_align/{sample}.dedup.bam",
        bed="results/03_peaks/consensus.bed",
    output:
        "results/qc/replicates/counts/{sample}.txt",
    params:
        pe=lambda wc: "-p --countReadPairs" if is_pe(wc.sample) else "",
    threads: threads_for("default")
    log:
        "<logs>/peak_counts/{sample}.log",
    conda:
        ENV_MACS
    shell:
        """
        awk -v OFS='\\t' 'BEGIN {{ print "GeneID", "Chr", "Start", "End", "Strand" }}
            {{ print $1 ":" $2 + 1 "-" $3, $1, $2 + 1, $3, "+" }}' {input.bed} > {output}.saf
        featureCounts -F SAF -a {output}.saf {params.pe} -T {threads} -o {output} {input.bam} > {log} 2>&1
        rm {output}.saf
        """


rule replicate_qc:
    """peak 重叠（Jaccard）、相关性、PCA、性别；结果为 tsv，绘图见 replicate_plots。"""
    input:
        counts=expand("results/qc/replicates/counts/{sample}.txt", sample=SAMPLES),
        peaks=[peak_file(s) for s in SAMPLES],
        reproducible=expand("results/03_peaks/reproducible/{group}.bed", group=REP_GROUPS),
        consensus="results/03_peaks/consensus.bed",
        sex=(
            expand("results/qc/replicates/sex/{sample}.tsv", sample=SAMPLES)
            if REPLICATES["sex_check"]
            else []
        ),
    output:
        samples="results/qc/replicates/samples.tsv",
        groups="results/qc/replicates/groups.tsv",
        overlap="results/qc/replicates/peak_overlap.tsv",
        venn="results/qc/replicates/peak_venn.tsv",
        correlation="results/qc/replicates/correlation.tsv",
        pca="results/qc/replicates/pca.tsv",
        counts="results/qc/replicates/peak_counts.tsv",
        info="results/qc/replicates/info.json",
    params:
        samples=SAMPLES,
        groups=GROUPS,
        rep_groups=REP_GROUPS,
        method=REPLICATES["method"],
        pca_top=REPLICATES["pca_top"],
        min_reps=REPLICATES["min_reps"],
        y_frac=REPLICATES["y_frac"],
        y_region=REF["y_region"],
        genome=config["genome"],
        mapq=config["filter"]["mapq"],
    log:
        "<logs>/replicate_qc.log",
    conda:
        ENV_MACS
    script:
        "../scripts/replicate_qc.py"


rule replicate_plots:
    """R（ggplot2）绘图：PNG 300 dpi、PDF、SVG（嵌入 QC 报告）。"""
    input:
        samples="results/qc/replicates/samples.tsv",
        groups="results/qc/replicates/groups.tsv",
        overlap="results/qc/replicates/peak_overlap.tsv",
        venn="results/qc/replicates/peak_venn.tsv",
        correlation="results/qc/replicates/correlation.tsv",
        pca="results/qc/replicates/pca.tsv",
        info="results/qc/replicates/info.json",
    output:
        replicate_figure_files(),
    params:
        figures=replicate_figures(),
        venn_groups=venn_groups(),
        outdir="results/qc/replicates/figures",
    log:
        "<logs>/replicate_plots.log",
    conda:
        ENV_MAIN
    script:
        "../scripts/replicate_plots.R"
