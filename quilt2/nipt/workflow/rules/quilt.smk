
rule quilt_prepare_regular:
    input:
        hap=rules.subset_refpanel_by_region.output.hap,
        leg=rules.subset_refpanel_by_region.output.leg,
        sam=rules.subset_refpanel_by_region.output.sam,
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "prep_regular",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_regular",
        nGen=config["quilt1"]["nGen"],
        buffer=config["quilt1"]["buffer"],
        lowram=config["quilt1"]["lowram"],
        impute_rare_common=config["quilt1"]["impute_rare_common"],
        rare_af_threshold=config["quilt1"]["rare_af_threshold"],
        gmap=if_use_quilt_map_in_refpanel,
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "prep_regular",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
        {params.time} -v QUILT_prepare_reference.R \
            --genetic_map_file='{params.gmap}' \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --reference_sample_file={input.sam} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --outputdir={params.outdir} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --reference_sample_file={input.sam} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --outputdir={params.outdir} \
        ; fi
        ) &> {log}
        """


rule quilt_run_regular:
    input:
        bam=rules.bamlist.output.bam,
        ff=rules.bamlist.output.ff,
        rdata=rules.quilt_prepare_regular.output,
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.regular.{chrom}.{start}.{end}.vcf.gz",
        ),
    params:
        time=config["time"],
        N="quilt_run_regular",
        method=config["quilt1"]["method"],
        nGen=config["quilt1"]["nGen"],
        buffer=config["quilt1"]["buffer"],
        Ksubset=config["quilt1"]["Ksubset"],
        nGibbsSamples=config["quilt1"]["nGibbsSamples"],
        n_seek_its=config["quilt1"]["n_seek_its"],
        block_gibbs=config["quilt1"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt1"]["small_ref_panel_gibbs_iterations"],
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.regular.{chrom}.{start}.{end}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ {params.method} == "nipt" ];then \
        {params.time} -v QUILT.R \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bam} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --fflist={input.ff} \
            --method="nipt" \
            --Ksubset={params.Ksubset} \
            --Knew={params.Ksubset} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --output_filename={output} \
        ; else \
        {params.time} -v QUILT.R \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bam} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --method="diploid" \
            --Ksubset={params.Ksubset} \
            --Knew={params.Ksubset} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --output_filename={output} \
        ; fi
        ) &> {log}
        """



rule quilt_ligate_regular:
    input:
        get_quilt_regular_output,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.regular.{chrom}.vcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT1,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.{ff}f.regular.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.regular.{chrom}.vcf.gz.llog",
        ),
    params:
        N="quilt_ligate_regular",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        ( \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type z --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ) &> {log}
        """

