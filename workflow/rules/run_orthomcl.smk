rule reformat_fasta:
    input:
        proteom = "resources/{specie}_proteins.fa"
    output:
        proteom_reformat = "results/OrthoMCL/reformat/{specie}_proteins_reformat.fa"
    log:
        "logs/reformat_fasta/{specie}.log"
    conda:
        get_conda("orthomcl")
    threads: 1
    shell:
        """
        awk -v OFS='\t' '{{if ($0 ~ /^>/) print $1; else print $0;}}' {input.proteom} | 
        awk -F'|' '{{if ($0 ~ /^>/) print \">\" $3; else print $0;}}' | 
        sed  's/|>/|/g' |sed 's/\*//g' > {output.proteom_reformat}
        """

rule clean_proteins:
    input:
        rules.reformat_fasta.output
    output:
        link = temp("results/OrthoMCL/compliantFasta/{specie}.temp") # permet juste de passer à la règle suivante, le vrai output c'est results/OrthoMCL/compliantFasta
    log:
        "logs/clean_proteins/{specie}.log"
    conda:
        get_conda("orthomcl")
    threads:1
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
        get_conda("orthomcl")
    log:
        "logs/filter_proteins.log"
    threads:1
    shell:
        """
        cd results/OrthoMCL/ ; 
        orthomclFilterFasta compliantFasta 10 20
        """

rule make_blast:
    input:
        good_proteins = rules.filter_proteins.output.good_proteins
    output:
        bank = "results/OrthoMCL/blastBank/blastBank.pdb"
    conda:
        get_conda("orthomcl")
    log:
        "logs/make_blast.log"
    params:
        bank = "results/OrthoMCL/blastBank/blastBank"
    threads:1
    shell:
        """
        makeblastdb -in {input.good_proteins} -dbtype prot -parse_seqids -out {params.bank}
        """

def get_part(number):
    len_number = len(str(number))
    number = int(number) + 1
    numbers = list(range(1, number))
    zero_filled_numbers = list()
    for i in numbers:
        i = str(i)
        zero_filled_numbers.append(i.zfill(len_number))
    
    return zero_filled_numbers

rule split_proteins:
    input:
        good_proteins = rules.filter_proteins.output.good_proteins
    output:
        split_proteins = expand("results/OrthoMCL/splited_proteins/goodProteins.part-{part}.fasta", part = get_part(config['split_proteins']))
    params:
        number_split = config['split_proteins']
    conda:
        get_conda("orthomcl")
    log:
         "logs/split_proteins.log"
    threads:20
    shell:
        """
        fasta-splitter --n-parts {params.number_split} {input.good_proteins} --out-dir results/OrthoMCL/splited_proteins
        """

rule blast:
    input:
        bank = rules.make_blast.output.bank,
        query = "results/OrthoMCL/splited_proteins/goodProteins.part-{part}.fasta"
    output:
        goodProteins = "results/OrthoMCL/blast/goodProteins_blast.part-{part}.tsv"    
    conda:
        get_conda("orthomcl")
    log:
        "logs/blast_part-{part}.log"
    params:
        bank = "results/OrthoMCL/blastBank/blastBank"
    threads:300
    shell:
        """
        blastp -db {params.bank} -query {input.query} -outfmt 6 -out {output.goodProteins}
        """


rule merge_and_convert:
    input:
        blast_split = expand("results/OrthoMCL/blast/goodProteins_blast.part-{part}.tsv", part = get_part(config['split_proteins']))
    output:
        blast_all = "results/OrthoMCL/blast/blast_results.tsv",
        blast_mysql = "results/OrthoMCL/blast/similarSequences.txt"
    conda:
        get_conda("orthomcl")
    log:
        "logs/merge_and_convert.log"
    threads:10
    shell:
        """

        orthomclBlastParser {output.blast_all} results/OrthoMCL/compliantFasta/ >> {output.blast_mysql}
        """

rule orthomcl_db:
    input:
        blast_mysql = "results/OrthoMCL/blast/similarSequences.txt"
    output:
        orthomcl_config = "results/OrthoMCL/orthomcl/orthomcl.config",
    conda:
        get_conda("orthomcl")
    threads:20   
    log:
        "logs/orthomcl_db.log"
    params:
        dbLogin = config['config_orthomcl']['dbLogin'],
        dbPassword = config['config_orthomcl']['dbPassword'],
        localHost = config['config_orthomcl']['localHost']
    shell:
        """
        CompEvo=$(pwd)
        cd results/OrthoMCL/

        echo 'dbVendor=mysql
dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1:{params.localHost}:3306
dbLogin={params.dbLogin}
dbPassword={params.dbPassword}
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE' > $CompEvo/{output.orthomcl_config}

        # orthomclInstallSchema $CompEvo/{output.orthomcl_config}
        # orthomclLoadBlast $CompEvo/{output.orthomcl_config} $CompEvo/{input.blast_mysql}
        """

rule compute_pairwise_relationships:
    input:
        orthomcl_config = rules.orthomcl_db.output.orthomcl_config
    output:
        orthologs = "results/OrthoMCL/orthomcl/pairs/orthologs.txt",
        inparalogs = "results/OrthoMCL/orthomcl/pairs/inparalogs.txt",
        coorthologs = "results/OrthoMCL/orthomcl/pairs/coorthologs.txt",
        mclInput = "results/OrthoMCL/orthomcl/mclInput"
    conda:
        get_conda("orthomcl")
    log:
        "logs/compute_pairwise_relationships.log"
    shell:
        """
        mkdir -p results/OrthoMCL/orthomcl/
        cd results/OrthoMCL/orthomcl/

        orthomclPairs orthomcl.config pairs.log cleanup=no
        orthomclDumpPairsFiles orthomcl.config
        """


rule clustering:
    input:
        # mclInput = ancient(rules.compute_pairwise_relationships.output.mclInput),
        # good_proteins = ancient(rules.filter_proteins.output.good_proteins)
        mclInput = "results/OrthoMCL/orthomcl/mclInput",
        good_proteins = "results/OrthoMCL/goodProteins.fasta"
    output:
        ortholog_groups = "results/OrthoMCL/clustering/inflation_{inflation}/groups_{inflation}.txt",
        named_groups = "results/OrthoMCL/clustering/inflation_{inflation}/named_groups_{inflation}.txt",
        named_groups_freq = "results/OrthoMCL/clustering/inflation_{inflation}/named_groups_{inflation}_freq.txt",
        sco_list = "results/OrthoMCL/clustering/inflation_{inflation}/scos_list_{inflation}.txt",
        sco_groups = "results/OrthoMCL/clustering/inflation_{inflation}/scos_{inflation}.ids",
        named_sco_groups = "results/OrthoMCL/clustering/inflation_{inflation}/named_groups_{inflation}_scos.txt",
        orthogroups_seq = directory("results/OrthoMCL/clustering/inflation_{inflation}/orthogroups_{inflation}")
    params:
        CopyNumberGen = "workflow/scripts/CopyNumberGen.sh",
        ExtractSCOs = "workflow/scripts/ExtractSCOs.sh",
        ExtractSeq = "workflow/scripts/ExtractSeq.sh"
    conda:
        get_conda("orthomcl")
    log:
        "logs/clustering_{inflation}.log"
    shell:
        """
        CompEvo=$(pwd)

        # MCL program (Dongen 2000) will be used to cluster the pairs extracted in the previous steps to determine ortholog groups.
        mcl {input.mclInput} --abc -I {wildcards.inflation} -o {output.ortholog_groups}

        # Name the groups called by mcl program.
        orthomclMclToGroups OG{wildcards.inflation}_ 1000 < {output.ortholog_groups} > {output.named_groups}

        # Frequency table for each ortholog
        {params.CopyNumberGen} {output.named_groups} > {output.named_groups_freq}

        # Select only 1:1 orthologs
        {params.ExtractSCOs} {output.named_groups_freq} > {output.sco_list}

        # Create a list of orthogroups that have SCOs
        cut -f 1 {output.sco_list} | grep -v "OG_name" | sed 's/$/:/' > {output.sco_groups}
        grep -Fw -f {output.sco_groups} {output.named_groups} > {output.named_sco_groups}

        # Extract the sequences
        {params.ExtractSeq} -o {output.orthogroups_seq} {output.named_sco_groups} {input.good_proteins}

        """


def inflation():
    inflation = list(np.arange(1,1.55,0.1)) + list(np.arange(2,6.5,0.5))
    for i in range(0,len(inflation)):
        inflation[i] = str(np.round(inflation[i],1))
    return inflation

def get_outgroup(INFOS):
    outgroup = ""
    for key, value in INFOS.items():
        if 'outgroup' in value:
            if value['outgroup'] == True:
                outgroup += value['code'] + '_'
    outgroup += '_'
    outgroup = outgroup.replace('__','')
    return outgroup


rule best_orthogroup:
    input:
        all_table_in = "results/OrthoMCL/clustering/inflation_{best_inflation}/named_groups_{best_inflation}_freq.txt",
        SCO_table_in = "results/OrthoMCL/clustering/inflation_{best_inflation}/scos_list_{best_inflation}.txt"
    output:
        all_table_out = "results/OrthoMCL/clustering/inflation_{best_inflation}/named_groups_{best_inflation}_freq_bestOG.tsv", 
        SCO_table_out = "results/OrthoMCL/clustering/inflation_{best_inflation}/scos_list_{best_inflation}_bestOG.tsv"
    conda:
        get_conda("orthomcl")
    params:
        findBestOutgroup = "workflow/scripts/findBestOutgroup.py",
        # drawOrthologyGroups = "workflow/scripts/drawOrthologyGroups.R",
        ougroups = get_outgroup(INFOS)
    log:
        "logs/best_orthogroup/best_orthogroup_{best_inflation}.log"
    shell: 
        """
        python {params.findBestOutgroup} --inputAllFreq {input.all_table_in} --outputAllFreq {output.all_table_out} \
            --inputScoFreq {input.SCO_table_in} --outputScoFreq {output.SCO_table_out} --outgroups {params.ougroups}
        """

rule infos_to_json:
    output:
        "resources/infos.json"
    params:
        infos = INFOS
    run:
        with open(output[0], "w") as outfile:
            json.dump(params[0], outfile)

# trier les protéines dans les fasta pour les mêmes dans le bon ordre
rule organize_fasta:
    input:
        # best_inflation = rules.best_inflation.output.best_inflation,
        fasta_nr = "results/OrthoMCL/clustering/inflation_{best_inflation}/orthogroups_{best_inflation}/OG{best_inflation}_{group}.fa", 
        json = rules.infos_to_json.output
    output:
        fasta_r = "results/OrthoMCL/clustering/inflation_{best_inflation}/orthogroups_{best_inflation}_reformat/OG{best_inflation}_{group}_reformat.fa"
    conda:
        get_conda("orthomcl")
    params:
        organizeFasta = "workflow/scripts/organizeFasta.py",
    log:
        "logs/best_inflation_{best_inflation}_{group}.log"
    shell: 
        """
        echo python {params.organizeFasta} -i {input.fasta_nr} -o {output.fasta_r} -j {input.json} 
        """

# rule maaft:
#     input:
#         fasta_r = rules.organize_fasta.output.fasta_r
#     output:
#         alignment = ""

# rule concatenate_alignments:
#     input:
#         rules.maaft.output.alignement
#     output:
#         ""

# rule raxml:
#     input:
#         rules.concatenate_alignments.output
#     output:

# rule draw_tree:
#     input:
#         rules.raxml.output
#     output:
