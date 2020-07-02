task tsne{
	input{
		File    anndata
		String  project_name
		Float   perplexity         = 30
		Float   early_exaggeration = 12
		Float   learning_rate      = 1000
		String?  memory ="4G"
		String?  cpu="2"
	}
	command <<<
		set -e
		set -o pipefail
		python << code
		import scanpy as sc
		infile   = "~{anndata}"
		adata = sc.read(infile)
		kwargs = {
			"perplexity"  : ~{perplexity},
			"early_exaggeration"    : ~{early_exaggeration},
			"learning_rate"     : ~{learning_rate}
			}
		sc.tl.tsne(adata, **kwargs)
		adata.obsm.to_df()[['X_tsne1','X_tsne2']].to_csv('~{project_name}_Scanpy_X_tsne.csv')
		adata.write("~{project_name}_tsne.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_tsne.h5ad"
  }
}
