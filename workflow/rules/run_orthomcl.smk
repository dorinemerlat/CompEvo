rule reformat_fasta:
    input:
        proteom = "resources/{specie}_proteins.fa"
    output:
        proteom_reformat = "results/OrthoMCL/reformat/{specie}_proteins_reformat.fa"
    log:
        "logs/reformat_fasta/{specie}.log"
    conda:
        get_conda("orthomcl")
        # "orthomcl"
    threads: 20
    shell:
        """
        awk -v OFS='\t' '{{if ($0 ~ /^>/) print $1; else print $0;}}' {input.proteom} | 
        awk -F'|' '{{if ($0 ~ /^>/) print \">\" $3; else print $0;}}' | 
        sed  's/|>/|/g' > {output.proteom_reformat}
        """

rule clean_proteins:
    input:
        rules.reformat_fasta.output
    output:
        link = temp("results/OrthoMCL/compliantFasta/{specie}.temp"), # permet juste de passer à la règle suivante, le vrai output c'est results/OrthoMCL/compliantFasta
    log:
        "logs/clean_proteins/{specie}.log"
    conda:
        get_conda("orthomcl")
    threads:20
    params:
        specie="{specie}",
        code_taxon = get_code_taxon
    shell:
        """
        if [ ! -d results/OrthoMCL/compliantFasta ] ; then
           mkdir results/OrthoMCL/compliantFasta
           echo '1'
        fi
        cd results/OrthoMCL/compliantFasta
        # taxon_code=$(echo {params.specie} |awk -F'-' '{{print toupper(substr($1,1,2) substr($2,1,2))}}') 
        echo {params.code_taxon}
        orthomclAdjustFasta {params.code_taxon} ../../../{input} 1
        touch ../../../{output.link}
        """

        
rule filter_proteins:
    input:
        link = expand("results/OrthoMCL/compliantFasta/{specie}.temp", specie=SPECIES)
    output:
        good_proteins = "results/OrthoMCL/goodProteins.fasta"
    conda:
        # get_conda("orthomcl")
        "orthomcl"
    log:
        "logs/filter_proteins/filter_proteins.log"
    threads:20
    shell:
        """
        cd results/OrthoMCL/ ; 
        orthomclFilterFasta compliantFasta 10 20
        """

rule make_blast:
    input:
        good_proteins = rules.filter_proteins.output
    output:
        bank = "results/OrthoMCL/bank_blast.pdb"
    conda:
        get_conda("orthomcl")
    log:
        "logs/make_blast/make_blast.log"
    params:
        bank = "results/OrthoMCL/bank_blast"
    threads:100
    shell:
        """
        makeblastdb -in {input.good_proteins} -dbtype prot -parse_seqids -out {params.bank}
        """
        
rule split_proteins:
    input:
        good_proteins = rules.filter_proteins.output
    output:
        split_proteins = "results/OrthoMCL/blast/goodProteins.part-{part}.fasta"
    conda:
        get_conda("orthomcl")
    log:
         "logs/split_proteins/split_proteins-{part}.log"
    threads:1
    shell:
        """
        fasta-splitter --n-parts 10000 {input.good_proteins} --out-dir results/OrthoMCL/blast
        """

rule blast:
    input:
        bank = rules.make_blast.output.bank,
        query = rules.split_proteins.output.split_proteins
    output:
        goodProteins = "results/OrthoMCL/blast/goodProteins.part-{part}.tsv"    
    conda:
        get_conda("orthomcl")
    log:
        "logs/blast/goodProteins.part-{part}.log"
    params:
        bank = "results/OrthoMCL/bank_blast"
    threads:16
    shell:
        """
        blastp -db {params.bank} -query {input.query} -outfmt 6 -out {output.goodProteins}
        """

rule merge_and_convert:
    input:
        blast_split = dynamic("results/OrthoMCL/blast/goodProteins.part-{part}.tsv")
    output:
        blast_all = "results/OrthoMCL/blast_results.tsv"#,
        # blast_mysql = "results/OrthoMCL/similarSequences.txt"
    conda:
        get_conda("orthomcl")
    log:
        "logs/merge_and_convert/merge_and_convert.log"
    threads:1
    shell:
        """
        cat {input.blast_split} >> {output.blast_all}
        """
        # # orthomclBlastParser {output.blast_all} results/OrthoMCL/compliantFasta/ >> {output.blast_mysql}
        # """

# rule orthomcl_db:
#     input:
#         blast_mysql = rules.merge_and_convert.output.blast_mysql
#     output:
#         orthomcl_config = "results/OrthoMCL/orthomcl.config",
#         mysql_dir = directory("results/OrthoMCL/mysql")
#     conda:
#         # get_conda("orthomcl")
#         "orthomcl"
#     threads:20   
#     log:
#         "logs/orthomcl_db/orthomcl_db.log"
#     params:
#         dbLogin = config['config_orthomcl']['dbLogin'],
#         dbPassword = config['config_orthomcl']['dbPassword'],
#         localHost = config['config_orthomcl']['localHost']
#     shell:
#         """
#         cd results/OrthoMCL/

#         echo "dbVendor=mysql
# dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1:localhost:{params.localHost}
# dbLogin={params.dbLogin}
# dbPassword={params.dbPassword}
# similarSequencesTable=SimilarSequences
# orthologTable=Ortholog
# inParalogTable=InParalog
# coOrthologTable=CoOrtholog
# interTaxonMatchView=InterTaxonMatch
# percentMatchCutoff=50
# evalueExponentCutoff=-5
# oracleIndexTblSpc=NONE
#         " > {output.orthomcl_config}

#         orthomclInstallSchema {output.orthomcl_config}

#         orthomclLoadBlast orthomcl.config {input.blast_mysql}
#         """

# rule compute_pairwise_relationships:
#     input:
#         rules.orthomcl_db.output.orthomcl_config
#     output:
#         orthologs = "results/OrthoMCL/pairs/potentialOrthologs.txt",
#         inparalogs = "results/OrthoMCL/pairs/potentialInparalogs.txt",
#         coorthologs = "results/OrthoMCL/pairs/potentialCoorthologs.txt"
#     conda:
#         # get_conda("orthomcl")
#         "orthomcl"
#     log:
#         "logs/compute_pairwise_relationships/compute_pairwise_relationships.log"
#     shell:
#         """
#         cd results/OrthoMCL
#         orthomclPairs orthomcl.config pairs.log cleanup=no
#         orthomclDumpPairsFiles orthomcl.config
#         """