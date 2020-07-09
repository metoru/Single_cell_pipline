version 1.0
workflow main{
	input{
	     File anndata
	     String project_name
	   }
	 call hvg{
	    input:
	       anndata=anndata,
	       project_name= project_name
	   }
	 output {
	  File        h5adfile=hvg.outputfile
	  Array[File] pngfile=hvg.pngfile
	  }
}	       
task hvg{
	input {	
		File     anndata
		Float    min_mean
		String	 project_name
		Float    max_mean
		Float    min_disp
		Float?   max_disp
		Float    max_value 
		Int?     n_top_genes
		String?  flavor 
		Boolean? subset
		String?  memory ="4G"
		String?  cpu="2"
	}
	command <<<
		set -e
		set -o pipefail
		python << code
		import scanpy as sc
		infile   = "~{anndata}"
		kwargs = {
			"min_mean" : ~{min_mean},
			"flavor"   : "~{default='seurat' flavor}", 
			"max_mean" : ~{max_mean},
			"min_disp" : ~{min_disp}
			}
		maxdisp=float('~{default='inf' max_disp}')
		adata = sc.read(infile)
		sc.pp.highly_variable_genes(adata,max_disp=maxdisp,subset=bool(~{true=1 false=0 subset}) ~{",n_top_genes=" +n_top_genes}, **kwargs)
		sc.pl.highly_variable_genes(adata, save=".png")
		adata = adata[:, adata.var.highly_variable]
		sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])# log1
		sc.pp.scale(adata, max_value=~{max_value})
		adata.write("~{project_name}_hvg.h5ad", compression="gzip") 
		code
	>>>
	runtime {
		memory: memory
		cpu: cpu
	}
  output {
    	File outputfile = "${project_name}_hvg.h5ad"
 	Array[File] pngfile    = glob("figures/*png")
  }
}

