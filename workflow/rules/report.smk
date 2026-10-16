# 汇总表、MultiQC


def qc_files():
    """汇总 QC：samtools stats（含片段长度分布）；TSS profile 由 config 控制。"""
    d = {"stats": expand("results/02_align/stats/{sample}.stats", sample=SAMPLES)}
    if config["qc"]["tss_profile"]:
        d["tss"] = "results/qc/deeptools/tss_profile.tsv"
    return d


rule summary:
    input:
        unpack(lambda wc: qc_files()),
        trim=[f"results/01_trim/reports/{s}_R1.fastq.gz_trimming_report.txt" for s in SAMPLES],
        bowtie2=expand("results/qc/bowtie2/{sample}.log", sample=SAMPLES),
        idxstats=expand("results/02_align/stats/{sample}.raw.idxstats", sample=SAMPLES),
        filtered=expand("results/02_align/stats/{sample}.filtered.flagstat", sample=SAMPLES),
        final=expand("results/02_align/stats/{sample}.flagstat", sample=SAMPLES),
        frip=expand("results/qc/frip/{sample}.tsv", sample=SAMPLES),
        peaks=[peak_file(s) for s in SAMPLES],
        contam=(
            expand("results/qc/contamination/{sample}.log", sample=SAMPLES)
            if CONTAM["enabled"]
            else []
        ),
    output:
        tsv="results/tables/summary.tsv",
        mqc="results/qc/summary_mqc.tsv",
    params:
        samples=SAMPLES,
        meta={
            s: {
                "layout": "PE" if is_pe(s) else "SE",
            }
            for s in SAMPLES
        },
        mito=MITO,
        tss_flank=config["qc"]["tss_flank"],
    log:
        "<logs>/summary.log",
    script:
        "../scripts/summary.py"


rule multiqc:
    input:
        unpack(lambda wc: qc_files()),
        fastqc=[
            f"results/qc/{d}/{s}_{r}{suffix}_fastqc.zip"
            for s in SAMPLES
            for r in reads_of(s)
            for d, suffix in (("fastqc_raw", ""), ("fastqc_trimmed", ".trimmed"))
        ],
        trim=[
            f"results/01_trim/reports/{s}_{r}.fastq.gz_trimming_report.txt"
            for s in SAMPLES
            for r in reads_of(s)
        ],
        bowtie2=expand("results/qc/bowtie2/{sample}.log", sample=SAMPLES),
        flagstat=expand("results/02_align/stats/{sample}.flagstat", sample=SAMPLES),
        peaks=[peak_file(s) for s in SAMPLES],
        summary="results/qc/summary_mqc.tsv",
        config=workflow.source_path("../assets/multiqc_config.yaml"),
    output:
        "results/qc/multiqc/multiqc_report.html",
    params:
        # peak 文件本身 MultiQC 不识别，传入的是 macs3 的 xls（不是规则输出，按路径传）
        files=lambda wc, input: [
            f
            for k, v in input.items()
            if k not in ("config", "peaks")
            for f in ([v] if isinstance(v, str) else v)
        ]
        + [f"results/03_peaks/macs3/{s}/{s}_peaks.xls" for s in SAMPLES],
        title="bulk ATAC-seq",
    log:
        "<logs>/multiqc.log",
    conda:
        ENV_MAIN
    shell:
        """
        multiqc -f -c {input.config} --title '{params.title}' -o results/qc/multiqc \
            -n multiqc_report.html {params.files} > {log} 2>&1
        """

