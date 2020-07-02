task marker{
	input{
		File    anndata
		String  project_name
		String? method  = 't-test'
		String?  groupby 
		Int?    n_genes = 100
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
			"n_genes"  : ~{n_genes},
			"method"   : "~{method}"
		}
		sc.tl.rank_genes_groups(adata,  "~{groupby}",**kwargs)
		sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save=".png")
		adata.write("~{project_name}_marker.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
 	Array[File] pngfile    = glob("figures/*png")
    	File outputfile = "${project_name}_marker.h5ad"
  }
}
