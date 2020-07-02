task pca{
	input {
		File     anndata
		Int?     n_comps
		String   svd_solver 
		Boolean  zero_center
		Boolean?  use_highly_variable
		String   project_name
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
		sc.pp.pca(adata, n_comps=~{n_comps},svd_solver="~{svd_solver}",zero_center=bool(~{true=1 false=0 zero_center}),\
			~{true=",use_highly_variable=True" false=",use_highly_variable=False" use_highly_variable}
			)
		sc.pl.pca_variance_ratio(adata, log=True, save=".png")
		adata.write("~{project_name}_pca.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_pca.h5ad"
 	Array[File] pngfile    = glob("figures/*png")
  }
}
