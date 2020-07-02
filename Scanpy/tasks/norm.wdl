task normalize{
	input {	
		File     anndata
		String	 project_name
		Float    target_sum
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
		sc.pp.normalize_total(adata, target_sum=~{target_sum})
		sc.pp.log1p(adata)
		adata.raw = adata
		adata.write("~{project_name}_norm.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_norm.h5ad"
  }
}

