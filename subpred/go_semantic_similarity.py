from rpy2.robjects import r, StrVector, packages, pandas2ri
import rpy2.robjects as ro


def get_semantic_similarities(
    go_terms, measure, organism, ont="MF", tcss_cutoff="NULL"
):
    # DOC: https://yulab-smu.top/biomedical-knowledge-mining-book/GOSemSim.html
    assert measure in {"Resnik", "Lin", "Rel", "Jiang", "TCSS", "Wang"}
    assert organism in {"yeast"}
    assert ont in {"CC", "MF", "BP"}

    packages.importr("GOSemSim")
    # add libraries for other organism here, also install via conda
    # alternative: use AnnotationForge to create custom org.db
    organism_to_orgdb = {"yeast": "org.Sc.sgd.db"}
    packages.importr(organism_to_orgdb[organism])
    packages.importr("Rcpp")
    r(
        """
        get_semantic_similarity_scores <- function(go_filter_input,measure, ont, orgDbName, tcss_cutoff=NULL){
            # also add libraries for other organism here:
            go_annot <- godata(
                OrgDb = orgDbName,
                ont=ont, 
                computeIC = TRUE, 
                processTCSS = TRUE,
                # keytype = "UNIPROT"
                cutoff=tcss_cutoff
            )
            
            # unique_go_terms <- unique(as.vector(go_annot@geneAnno[["GO"]]))
            # go_filter <- c()
            # for (go_term in go_filter_input){
            #     if (go_term %in% unique_go_terms){
            #         go_filter <- c(go_filter, go_term)
            #     }
            # }
            go_filter_input <- as.character(go_filter_input)
            
            values <- mgoSim(
                go_filter_input,
                go_filter_input,
                semData=go_annot,
                measure=measure,
                combine=NULL
            )
            
            return (as.data.frame(values))
        }
    """
    )

    matr = r["get_semantic_similarity_scores"](
        go_terms,
        measure=measure,
        ont=ont,
        orgDbName=organism_to_orgdb[organism],
        tcss_cutoff=ro.NULL if tcss_cutoff == "NULL" else tcss_cutoff,
    )

    with (ro.default_converter + pandas2ri.converter).context():
        df_go_similarity = ro.conversion.get_conversion().rpy2py(matr)
    df_go_similarity = df_go_similarity.loc[sorted(go_terms), sorted(go_terms)]
    return df_go_similarity


# # https://rpy2.github.io/doc/latest/html/introduction.html
# # https://berkeley-scf.github.io/tutorial-parallelization/parallel-R.html

# from tempfile import NamedTemporaryFile
# from rpy2.robjects import r, StrVector, packages, pandas2ri
# import rpy2.robjects as ro
# import numpy as np
# import pandas as pd
# from subpred.util import load_df

# def get_go_data_df(
#     organism_ids: set,
#     aspect: str,
# ):
#     """
#         The GoSemSim package expects data from AnnotationDBI.
#         However, this means that we can only calculate similarities between GO terms from the same organism, not on the entire GO.
#         We fix this by creating this object directly, from our go annotation pickle, analogous to the godata function in GoSemSim_summary_.
#         The original code does not use the GO ancestors (GOALL), only direct annotations (GO).
#         It adds the ancestors later during calculation if they are needed.
#         If an AnnotationDbi object is available then we use that.
#     Args:
#         organism_ids (set): Organism ids to take annotated proteins from
#         aspect (str): Sub-ontology, one of C, F, P. Defaults to F

#     Returns:
#         pd.DataFrame: Annotation dataframe compatible with get_go_semantic_similarities
#     """
#     df_uniprot_org = load_df("uniprot")
#     df_uniprot_org = df_uniprot_org[
#         df_uniprot_org.organism_id.isin(organism_ids)
#     ]  # TODO filter for more?
#     protein_set = set(df_uniprot_org.index)

#     df_go_data_semsim = load_df("go")
#     df_go_data_semsim = df_go_data_semsim[
#         (df_go_data_semsim.Uniprot.isin(protein_set))
#         # & (~df_go_data_semsim.evidence_code.isin(evidence_codes_exclude))
#         & (df_go_data_semsim.aspect == aspect)
#     ].reset_index(drop=True)
#     aspect_to_orgdb_aspect = {"F": "MF", "C": "CC", "P": "BP"}
#     df_go_data_semsim["aspect"] = df_go_data_semsim.aspect.map(aspect_to_orgdb_aspect)
#     df_go_data_semsim = df_go_data_semsim.drop("qualifier", axis=1).rename(
#         columns={
#             "Uniprot": "UNIPROT",
#             "go_id": "GO",
#             "evidence_code": "EVIDENCE",
#             "aspect": "ONTOLOGY",
#         }
#     )
#     return df_go_data_semsim


# def get_go_semantic_similarities(
#     organism_ids: set,
#     aspect: str,
#     go_id_subset: np.array,
#     method="Wang",
#     conversion_method: str = "internal",
#     parallelization_method:str="sequential"
# ):
#     """Wrapper for the GOSemSim R package, with support for efficient parallel processing

#     Args:
#         organism_ids (set): Organism ids to take annotated proteins from
#         aspect (str): Sub-ontology, one of C, F, P
#         go_id_subset (np.array): go ids for which to calculate pairwise similarities
#         method (str, optional):
#             Algorithm for calculation similarity.
#             Options: "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang"
#             Defaults to "Wang". Other methods require IC computation, only implemented for yeast
#         parallelization_method (str, optional):
#             If "multisession", clones the R session for every thread.
#                 Slower start and much higher memory consumption. Less likely to crash
#             If "multicore", forks R process and allows parallel access to the dataframes.
#                 Global vars *are* copied if they are modified by worker, but we only read them anyways.
#                 Much faster and lower memory consumption, not available on Windows.
#                 Can cause crashes for some values of parameter "method"
#             If "sequential": No parallelization
#             Defaults to "sequential".
#         conversion_method (str, optional):
#             If "temp_file", write results to tmp file in R and then read that file with pandas, deleting the file afterwards.
#             If "internal", uses internal conversion methods between rpy2 and pandas.
#             Which method is faster depends on the machine and the dataset size.
#             Defaults to "internal".
#     Raises:
#         ValueError: Parameter not recognized

#     Returns:
#         _type_: _description_
#     """

#     packages.importr("GOSemSim")
#     packages.importr("future.apply")
#     packages.importr("Rcpp")

#     r(
#         """
#         read_go_data <- function(path){
#             goa <- read.table(path)
#             kk <- unique(goa$UNIPROT)
#             d <- new("GOSemSimDATA", keys = kk, ont = "MF", geneAnno = goa)
#             return(d)
#         }
#         similarities <- function(go1, go2){
#             goSim(go1, go2, semData=d, measure="Wang")
#         }
#         similarities_vec <- function(go_term, go_vec, measure, go_data){
#             similarities_list <- lapply(go_vec, goSim, GOID2=go_term, semData=go_data, measure=measure)
#             return(unlist(similarities_list))
#         }

#         similarities_par <- function(go_vec, method, go_data, plan="sequential"){
#             if(plan == "sequential"){
#                 plan(sequential)
#             } else if (plan == "multicore"){
#                 plan(multicore)
#             } else if (plan == "multisession"){
#                 plan(multisession)
#             }
#             results <- future_lapply(go_vec, similarities_vec, go_vec=go_vec, measure=method, go_data=go_data)
#             # results <- lapply(go_vec, similarities_vec, go_vec=go_vec, measure=method, go_data=go_data)
#             results <- as.data.frame(results, row.names = go_vec, col.names = go_vec)
#             return(results)
#         }
#         write_tmp <- function(df, file_name){
#             write.table(df, file_name, sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
#         }
#     """
#     )
#     if organism_ids == {559292}:  # s. cerevisiae
#         packages.importr("org.Sc.sgd.db")
#         godata_func = r["godata"]
#         go_data_r = godata_func(
#             "org.Sc.sgd.db",
#             # keytype="ENTREZID",
#             ont="MF",
#             computeIC=method != "Wang",
#             processTCSS=method == "TCSS",
#         )
#     else:
#         # Workaround for other datasets:Transform the data to what the R package expects
#         print("Warning: organism not found. creating custom GO dataset...")
#         df_go_data_semsim = get_go_data_df(organism_ids=organism_ids, aspect=aspect)

#         tmp_file = NamedTemporaryFile(suffix=".tsv")
#         df_go_data_semsim.to_csv(tmp_file, sep="\t")
#         df_goa_path = tmp_file.name
#         go_data_r = r["read_go_data"](df_goa_path)

#         tmp_file.close()

#     goSimVec = r["similarities_par"]
#     vector = (
#         np.unique(go_id_subset)
#         if go_id_subset is not None
#         else df_go_data_semsim.GO.unique()
#     )
#     vector_r = StrVector(vector)
#     results_r = goSimVec(vector_r, method, go_data_r, plan=parallelization_method)

#     if conversion_method == "internal":
#         # df_go_similarity = ro.conversion.rpy2py(results_r)
#         with (ro.default_converter + pandas2ri.converter).context():
#             df_go_similarity = ro.conversion.get_conversion().rpy2py(results_r)
#     elif conversion_method == "temp_file":
#         # write tmp file then read it with pandas. converting the R object to python took longer than the calculations
#         write_table = r["write_tmp"]
#         tmp_file = NamedTemporaryFile()
#         write_table(results_r, tmp_file.name)
#         df_go_similarity = pd.read_table(tmp_file)
#         tmp_file.close()
#     else:
#         raise ValueError("Invalid dataframe conversion method.")

#     df_go_similarity.columns = df_go_similarity.columns.to_series().str.replace(
#         ".", ":"
#     )
#     return df_go_similarity
