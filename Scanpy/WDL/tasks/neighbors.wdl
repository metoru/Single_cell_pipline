task neighbors{
	input {
		File      anndata
		Int       n_neighbors
		Int?      n_pcs
		String?   use_rep
		Boolean   knn
		String    method
		String?    key_added
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
		kwargs = {
			"n_neighbors"  : ~{n_neighbors},
			"method"       : "~{method}"
			}
		sc.pp.neighbors(adata, knn=bool(~{true=1 false=0 knn}) 
			~{   ",n_pcs="  + n_pcs} \
			~{ ",use_rep="  + use_rep}\
			~{",key_added="  + key_added}\
			,**kwargs)
		adata.write("~{project_name}_neighbor.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_neighbor.h5ad"
  }
}

