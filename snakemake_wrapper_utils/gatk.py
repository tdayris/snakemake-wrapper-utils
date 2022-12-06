import sys


def get_gatk4_opts(
    snakemake,
    parse_ref=True,
    parse_ref_dict=True,
    parse_bam_index=True,
    parse_bam_md5_digest=True,
    parse_tbi=True,
    parse_vcf_md5_digest=True,
    parse_intervals=True,
    parse_argument_file=True,
):
    """Obtain gatk4 opts from input, output, params and resources"""
    gatk_opts = ""
    extra = snakemake.params.get_extra("extra", "")

    ######################
    ### Reference file ###
    ######################
    if parse_ref:
        if "-R" in extra or "--reference" in extra:
            sys.exit(
                "You have specified reference path (`-R/--reference`) in params.extra; this is automatically infered from `ref` input file"
            )

        reference = snakemake.input.get("ref")
        if reference:
            gatk_opts += f" --reference {reference}"

    ###########################
    ### Sequence dictionary ###
    ###########################
    if parse_ref_dict:
        if "--sequence-dictionary" in extra:
            sys.exit(
                "You have specified reference dictionary path (`--sequence-dictionary`) in params.extra; this is automatically infered from `ref_dict` input file"
            )

        ref_dict = snakemake.output.get("ref_dict")
        if bam_index:
            gatk_opts += f" --sequence-dictionary {ref_dict}"

    #######################
    ### Write bam index ###
    #######################
    if parse_bam_index:
        if "-OBI" in extra or "--create-output-bam-index" in extra:
            sys.exit(
                "You have specified output bam index creation (`-OBI/--create-output-bam-index`) in params.extra; this is automatically infered from `bai` output file"
            )

        bam_index = snakemake.output.get("bai")
        if bam_index:
            gatk_opts += " --create-output-bam-index true"

    ########################
    ### Write md5 digest ###
    ########################
    if parse_bam_md5_digest:
        if "-OBM" in extra or "--create-output-bam-md5" in extra:
            sys.exit(
                "You have specified output bam MD5 digest creation (`-OBM/--create-output-bam-md5`) in params.extra; this is automatically infered from `bam_md5` output file"
            )

        bam_md5 = snakemake.output.get("bam_md5")
        if bam_md5:
            gatk_opts += " --create-output-bam-md5 true"

    ########################
    ### Write vcf index  ###
    ########################
    if parse_tbi:
        if "-OVI" in extra or "--create-output-variant-index" in extra:
            sys.exit(
                "You have specified output VCF index creation (`-OVI/--create-output-variant-index`) in params.extra; this is automatically infered from `tbi` output file"
            )

        bam_md5 = snakemake.output.get("tbi")
        if bam_md5:
            gatk_opts += " --create-output-variant-index true"

    #############################
    ### Write vcf MD5 digest  ###
    #############################
    if parse_vcf_md5_digest:
        if "-OVM" in extra or "--create-output-variant-md5" in extra:
            sys.exit(
                "You have specified output VCF MD5 digest creation (`-OVM/--create-output-variant-md5`) in params.extra; this is automatically infered from `vcf_md5` output file"
            )

        bam_md5 = snakemake.output.get("vcf_md5")
        if bam_md5:
            gatk_opts += " --create-output-variant-md5 true"

    #########################
    ### Use interval file ###
    #########################
    if parse_intervals:
        if "--intervals" in extra or "-L" in extra:
            sys.exit(
                "You have specified input interval file (`-L/--intervals`) in params.extra; this is automatically infered from `intervals` input file"
            )

        intervals = snakemake.input.get("intervals")
        if intervals:
            gatk_opts += f" --intervals {intervals}"

    #########################
    ### Use argument file ###
    #########################
    if parse_argument_file:
        if "--arguments_file" in extra:
            sys.exit(
                "You have specified input argument file (`--arguments_file`) in params.extra; this is automatically infered from `arguments` input file"
            )

        args = snakemake.input.get("arguments")
        if intervals:
            gatk_opts += f" --intervals {args}"

    return gatk_opts
