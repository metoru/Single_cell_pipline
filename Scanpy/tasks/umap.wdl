task umap{
	input{
		File    anndata
		Float   min_dist
		Float   spread
		Float   alpha
		Float   gamma
		Int     n_components
		Int     negative_sample_rate
		String  init_pos
		String  project_name
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
			"min_dist"  : ~{min_dist},
			"spread"    : ~{spread},
			"alpha"     : ~{alpha},
			"gamma"     : ~{gamma},
			"n_components" : ~{n_components},
			"negative_sample_rate" : ~{negative_sample_rate},
			"init_pos"  : "~{init_pos}"
			}
		sc.tl.umap(adata, **kwargs)
		adata.obsm.to_df()[['X_umap1','X_umap2']].to_csv('~{project_name}_Scanpy_X_umap.csv')
		adata.write("~{project_name}_umap.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_umap.h5ad"
	File clustfile  = "${project_name}_Scanpy_X_umap.csv"
  }
}

