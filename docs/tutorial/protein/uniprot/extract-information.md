# From text


The human proteome is downloaded from [UniProt - Reviewed Swiss-Prot](https://www.uniprot.org/uniprotkb?query=Human&facets=reviewed%3Atrue%2Cmodel_organism%3A9606). As of Jule 27th, 2024, it has compiled 20435 expert-reviewed human genes that encode proteins. As UniProt sequences in [Text](https://www.uniprot.org/help/retrieve_sets#:~:text=UniProtKB%20entries%20are%20available%20in,%2Dcalled%20'canonical'%20sequence.) format contains information about the sequences, structures,and functions of proteins, users can extract as comprehensive information as possible from the Text file.

The following examples show how we use `pypropel` to achieve this functionality.

``` py linenums="1"
import pypropel as pp

df = pp.uniprot.from_text(
    text_fpn=to('data/uniprot/text/uniprotkb_Human_AND_reviewed_true_AND_m_2023_11_29.txt'),
    sv_json_fpn=to('data/uniprot/text/human.json'),
    sv_df_fpn=to('data/uniprot/text/human.txt'),
)
```

The output of the `pp.uniprot.from_text` function.
``` shell
ID	AC	DE	GN	Ensembl_G_id	Ensembl_T_id	SQ	tms_nums	binding_nums	binding_sites	pdb_nums	pdb_ids	pdb_chains	pdb_rezs
CP2D7_HUMAN	A0A087X1C5	RecName: Full=Putative cytochrome P450 2D7 {ECO:0000305}; EC=1.14.14.1 {ECO:0000269|PubMed:15051713};	CYP2D7	NaN		MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQLRRRFGDVFSLQLAWTPVVVLNGLAAVREAMVTRGEDTADRPPAPIYQVLGFGPRSQGVILSRYGPAWREQRRFSVSTLRNLGLGKKSLEQWVTEEAACLCAAFADQAGRPFRPNGLLDKAVSNVIASLTCGRRFEYDDPRFLRLLDLAQEGLKEESGFLREVLNAVPVLPHIPALAGKVLRFQKAFLTQLDELLTEHRMTWDPAQPPRDLTEAFLAKKEKAKGSPESSFNDENLRIVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQRGRRVSPGCPIVGTHVCPVRVQQEIDDVIGQVRRPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSRDIEVQGFRIPKGTTLITNLSSVLKDEAVWKKPFRFHPEHFLDAQGHFVKPEAFLPFSAGRRACLGEPLARMELFLFFTSLLQHFSFSVAAGQPRPSHSRVVSFLVTPSPYELCAVPR	2	1	461	0			
PIOS1_HUMAN	A0A0B4J2F0	RecName: Full=Protein PIGBOS1 {ECO:0000305}; AltName: Full=PIGB opposite strand protein 1 {ECO:0000312|HGNC:HGNC:50696};	PIGBOS1	ENSG00000225973	ENST00000436697.3;ENST00000567948.1	MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS	1	0		0	
E2F8_HUMAN	A0AVK6	RecName: Full=Transcription factor E2F8; Short=E2F-8;	E2F8	ENSG00000129173	ENST00000250024.9;ENST00000527884.5;ENST00000620009.4	MENEKENLFCEPHKRGLMKTPLKESTTANIVLAEIQPDFGPLTTPTKPKEGSQGEPWTPTANLKMLISAVSPEIRNRDQKRGLFDNRSGLPEAKDCIHEHLSGDEFEKSQPSRKEKSLGLLCHKFLARYPNYPNPAVNNDICLDEVAEELNVERRRIYDIVNVLESLHMVSRLAKNRYTWHGRHNLNKTLGTLKSIGEENKYAEQIMMIKKKEYEQEFDFIKSYSIEDHIIKSNTGPNGHPDMCFVELPGVEFRAASVNSRKDKSLRVMSQKFVMLFLVSTPQIVSLEVAAKILIGEDHVEDLDKSKFKTKIRRLYDIANVLSSLDLIKKVHVTEERGRKPAFKWTGPEISPNTSGSSPVIHFTPSDLEVRRSSKENCAKNLFSTRGKPNFTRHPSLIKLVKSIESDRRKINSAPSSPIKTNKAESSQNSAPFPSKMAQLAAICKMQLEEQSSESRQKVKVQLARSGPCKPVAPLDPPVNAEMELTAPSLIQPLGMVPLIPSPLSSAVPLILPQAPSGPSYAIYLQPTQAHQSVTPPQGLSPTVCTTHSSKATGSKDSTDATTEKAANDTSKASASTRPGSLLPAPERQGAKSRTREPAGERGSKRASMLEDSGSKKKFKEDLKGLENVSATLFPSGYLIPLTQCSSLGAESILSGKENSSALSPNHRIYSSPIAGVIPVTSSELTAVNFPSFHVTPLKLMVSPTSVAAVPVGNSPALASSHPVPIQNPSSAIVNFTLQHLGLISPNVQLSASPGSGIVPVSPRIESVNVAPENAGTQQGRATNYDSPVPGQSQPNGQSVAVTGAQQPVPVTPKGSQLVAESFFRTPGGPTKPTSSSCMDFEGANKTSLGTLFVPQRKLEVSTEDVH	0	0		1	4YO2	A=110-341	3.07
...
```

The information extracted is shown as follows. 
``` json
{
    'ID': Uniprot ,
    'AC': Uniprot accession code,
    'DE': Description,
    'GN': gene names,
    'Ensembl_G_id': ensembl gene ID,
    'Ensembl_T_id': ensembl transcript ID,
    'SQ': sequence,
    'tms_nums': number of transmembrane segments,
    'binding_nums': numbers of binding sites,
    'binding_sites': PDB IDs of binding sites,
    'pdb_nums': number of PDB structures,
    'pdb_ids': IDS of PDB structures,
    'pdb_chains': IDs of PDB chains,
    'pdb_rezs': resolutions of PDB structures,
}
```

This function will also output a file in `json` format. It gives specific PDB IDs of the extracted transmembrane segments in `tmh_lower` and `tmh_upper`.
``` shell
{"CP2D7_HUMAN": {"pac": "A0A087X1C5", "gn": "CYP2D7", "ensembl_id": "NaN", "binding_num": 1, "pdb_num": 0, "cyto_lower": [24], "cyto_upper": [301], "extra_lower": [1, 323], "extra_upper": [2, 515], "tmh_lower": [3, 302], "tmh_upper": [23, 322], "nontmh_lower": [1, 24, 323], "nontmh_upper": [2, 301, 515], "seq": "MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQLRRRFGDVFSLQLAWTPVVVLNGLAAVREAMVTRGEDTADRPPAPIYQVLGFGPRSQGVILSRYGPAWREQRRFSVSTLRNLGLGKKSLEQWVTEEAACLCAAFADQAGRPFRPNGLLDKAVSNVIASLTCGRRFEYDDPRFLRLLDLAQEGLKEESGFLREVLNAVPVLPHIPALAGKVLRFQKAFLTQLDELLTEHRMTWDPAQPPRDLTEAFLAKKEKAKGSPESSFNDENLRIVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQRGRRVSPGCPIVGTHVCPVRVQQEIDDVIGQVRRPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSRDIEVQGFRIPKGTTLITNLSSVLKDEAVWKKPFRFHPEHFLDAQGHFVKPEAFLPFSAGRRACLGEPLARMELFLFFTSLLQHFSFSVAAGQPRPSHSRVVSFLVTPSPYELCAVPR"}, "PIOS1_HUMAN": {"pac": "A0A0B4J2F0", "gn": "PIGBOS1", "ensembl_id": "ENSG00000225973", "binding_num": 0, "pdb_num": 0, "cyto_lower": [26], "cyto_upper": [54], "extra_lower": [], "extra_upper": [], "tmh_lower": [5], "tmh_upper": [25], "nontmh_lower": [1, 26], "nontmh_upper": [4, 54], "seq": "MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS"},
...
}
```

!!! tip

    The Json file also returns cytoplasmic and extracellular segments of the extracted transmembrane segments, as indicated by `cyto_lower`, `cyto_upper`,  `extra_lower`, and `extra_upper`.  Their combined segments are non-transmembrane segments (`nontmh_lower` amd `nontmh_upper`).