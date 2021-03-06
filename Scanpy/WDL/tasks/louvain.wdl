task leiden{
	input{
		File    anndata
		String  project_name
		Float   resolution          = 1
		Boolean directed            = true
		Boolean   use_weights         = true
		Float   n_iterations       = -1
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
			"resolution"  : ~{resolution},
			"directed"    : bool(~{true=1 false=0 directed}),
			"use_weights"     : bool(~{true=1 false=0 use_weights}),
			"n_iterations"    : ~{n_iterations},
			}
		sc.tl.leiden(adata, **kwargs)
		if 'leiden' not in adata.obs.keys():
			raise KeyError('leiden is not a valid key')
		#adata.obs['leiden'].reset_index(level=0).rename(columns={'index': 'cells'}).to_csv('clust.txt', sep='\t', header=True, index=False)
		adata.obs[['n_genes','total_counts','leiden']].to_csv('~{project_name}_clust.txt', sep='\t', header=True, index=False)
		adata.write("~{project_name}_leiden.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_leiden.h5ad"
	File clustfile  =  "${project_name}_clust.txt"
  }
}


