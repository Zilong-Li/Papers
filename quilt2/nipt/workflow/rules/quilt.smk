
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


rule quilt_prepare_mspbwt:
    input:
        vcf=rules.subset_refpanel_by_region.output.vcf,
    output:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "prep_mspbwt",
            "RData",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_mspbwt",
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
        nGen=config["quilt2"]["nGen"],
        buffer=config["quilt2"]["buffer"],
        gmap=if_use_quilt_map_in_refpanel,
        lowram=config["quilt2"]["lowram"],
        impute_rare_common=config["quilt2"]["impute_rare_common"],
        rare_af_threshold=config["quilt2"]["rare_af_threshold"],
        nindices=config["quilt2"]["mspbwt-nindices"],
    log:
        lambda wildcards, output: output[0] + ".llog",
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
        {params.time} -v QUILT_prepare_reference.R \
            --genetic_map_file='{params.gmap}' \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --use_hapMatcherR={params.lowram} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_mspbwt=TRUE \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --mspbwt_nindices={params.nindices} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --use_hapMatcherR={params.lowram} \
            --nGen={params.nGen} \
            --use_mspbwt=TRUE \
            --rare_af_threshold={params.rare_af_threshold} \
            --impute_rare_common={params.impute_rare_common} \
            --mspbwt_nindices={params.nindices} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; fi \
        ) &> {log}
        """

rule quilt_run_mspbwt:
    input:
        vcf=rules.subset_refpanel_by_region.output.vcf,
        bam=rules.bamlist.output.bam,
        ff=rules.bamlist.output.ff,
        rdata=rules.quilt_prepare_regular.output,
    output:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.mspbwt.{chrom}.{start}.{end}.vcf.gz",
        ),
    params:
        time=config["time"],
        N="quilt_run_mspbwt",
        nGen=config["quilt2"]["nGen"],
        buffer=config["quilt2"]["buffer"],
        Ksubset=config["quilt2"]["Ksubset"],
        nGibbsSamples=config["quilt2"]["nGibbsSamples"],
        n_seek_its=config["quilt2"]["n_seek_its"],
        lowram=config["quilt2"]["lowram"],
        rare_af_threshold=config["quilt2"]["rare_af_threshold"],
        impute_rare_common=config["quilt2"]["impute_rare_common"],
        block_gibbs=config["quilt2"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt2"]["small_ref_panel_gibbs_iterations"],
        mspbwtM=config["quilt2"]["mspbwtM"],
        mspbwtL=config["quilt2"]["mspbwtL"],
    log:
        lambda wildcards, output: output[0] + ".llog",
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ {params.method} == "nipt" ];then \
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bam} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --fflist={input.ff} \
            --method="nipt" \
            --Knew={params.Ksubset} \
            --use_hapMatcherR={params.lowram} \
            --impute_rare_common={params.impute_rare_common} \
            --zilong=FALSE \
            --use_mspbwt=TRUE \
            --mspbwtM={params.mspbwtM} \
            --mspbwtL={params.mspbwtL} \
            --rare_af_threshold={params.rare_af_threshold} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --output_filename={output} \
        ; else \
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bam} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --method="diploid" \
            --Knew={params.Ksubset} \
            --use_hapMatcherR={params.lowram} \
            --impute_rare_common={params.impute_rare_common} \
            --zilong=FALSE \
            --use_mspbwt=TRUE \
            --mspbwtM={params.mspbwtM} \
            --mspbwtL={params.mspbwtL} \
            --rare_af_threshold={params.rare_af_threshold} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --output_filename={output} \
        ; fi
        ) &> {log}
        """
        
rule quilt_ligate_mspbwt:
    input:
        get_quilt_mspbwt_output,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.{ff}f.mspbwt.{chrom}.vcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT1,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.{ff}f.mspbwt.{chrom}.vcf.list",
            )
        ),
    log:
        lambda wildcards, output: output[0] + ".llog",
    params:
        N="quilt_ligate_mspbwt",
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
